using System;
using TMPro;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEngine.UI;

public class InflationAnimationModule : MonoBehaviour, IAnimationModule
{
    public string Id => "inflation";
    public string Title => "インフレーション";

    private const float BoxW = 4.2f;
    private const float BoxH = 10.5f;
    private const float MinLog = -31f;
    private const float MaxLog = 0f;
    private const float StepInterval = 0.04f;

    private static readonly Color[] Palette =
    {
        HexColor(0x636efa), HexColor(0xEF553B), HexColor(0x00cc96), HexColor(0xab63fa), HexColor(0xFFA15A),
        HexColor(0x19d3f3), HexColor(0xFF6692), HexColor(0xB6E880), HexColor(0xFF97FF), HexColor(0xFECB52)
    };

    private int _m;
    private int _t;
    private float _floor;
    private float _rmax;
    private float[] _wr;
    private float[] _wi;
    private float[] _fz;
    private bool _dataLoaded;

    private GameObject _sceneRoot;
    private Mesh _fiberMesh;
    private Mesh _markerMesh;
    private Vector3[] _markerVerts;
    private readonly Transform[] _billboards = new Transform[12];
    private int _billboardCount;

    private RectTransform _uiContainer;
    private Slider _slider;
    private TextMeshProUGUI _statusLabel;
    private TextMeshProUGUI _playBtnLabel;
    private InflationCurveGraphic _curveGraphic;

    private Camera _camera;
    private TMP_FontAsset _font;
    private Material _vertexColorMat;

    private int _currentStep;
    private bool _isPlaying;
    private float _playTimer;
    private bool _active;
    private bool _suppressSliderCallback;

    public void Activate(RectTransform panelRoot, Camera mainCamera, TMP_FontAsset font)
    {
        _camera = mainCamera;
        _font = font;
        _active = true;

        EnsureDataLoaded();
        if (!_dataLoaded) return;

        if (_camera != null)
        {
            _camera.clearFlags = CameraClearFlags.SolidColor;
            _camera.backgroundColor = Color.white;
            var urpData = _camera.GetUniversalAdditionalCameraData();
            if (urpData != null)
            {
                urpData.renderPostProcessing = false;
            }
            var orbit = _camera.GetComponent<OrbitCameraController>();
            if (orbit == null) orbit = _camera.gameObject.AddComponent<OrbitCameraController>();
            orbit.minDistance = 4f;
            orbit.maxDistance = 40f;
            orbit.SetView(38f, 14f, 18.0f, new Vector3(0f, -1.15f, 0f));
        }

        if (_sceneRoot == null)
        {
            BuildScene3D();
        }
        _sceneRoot.SetActive(true);

        if (_uiContainer == null)
        {
            BuildUI(panelRoot);
        }
        _uiContainer.gameObject.SetActive(true);

        SetStep(0);
        SetPlaying(true);
    }

    public void Deactivate()
    {
        _active = false;
        SetPlaying(false);
        if (_sceneRoot != null) _sceneRoot.SetActive(false);
        if (_uiContainer != null) _uiContainer.gameObject.SetActive(false);
    }

    private void Update()
    {
        if (!_active || !_dataLoaded) return;

        if (_camera != null)
        {
            Quaternion camRot = _camera.transform.rotation;
            for (int i = 0; i < _billboardCount; i++)
            {
                if (_billboards[i] != null) _billboards[i].rotation = camRot;
            }
            UpdateMarkersBillboard();
        }

        if (_isPlaying)
        {
            _playTimer += Time.deltaTime;
            while (_playTimer >= StepInterval)
            {
                _playTimer -= StepInterval;
                if (_currentStep + 1 < _t)
                {
                    SetStep(_currentStep + 1);
                }
                else
                {
                    SetPlaying(false);
                    break;
                }
            }
        }
    }

    private void EnsureDataLoaded()
    {
        if (_dataLoaded) return;

        TextAsset asset = Resources.Load<TextAsset>("InflationN40");
        if (asset == null)
        {
            Debug.LogError("[InflationAnimationModule] Resources/InflationN40.bytes not found.");
            return;
        }

        byte[] bytes = asset.bytes;
        int offset = 0;
        _m = BitConverter.ToInt32(bytes, offset); offset += 4;
        _t = BitConverter.ToInt32(bytes, offset); offset += 4;
        _floor = BitConverter.ToSingle(bytes, offset); offset += 4;

        int waveCount = _m * _t;
        _wr = new float[waveCount];
        _wi = new float[waveCount];
        _fz = new float[_t];

        Buffer.BlockCopy(bytes, offset, _wr, 0, waveCount * 4); offset += waveCount * 4;
        Buffer.BlockCopy(bytes, offset, _wi, 0, waveCount * 4); offset += waveCount * 4;
        Buffer.BlockCopy(bytes, offset, _fz, 0, _t * 4);

        float maxAbs = 0f;
        for (int i = 0; i < waveCount; i++)
        {
            float r = Mathf.Sqrt(_wr[i] * _wr[i] + _wi[i] * _wi[i]);
            if (r > maxAbs) maxAbs = r;
        }
        _rmax = Mathf.Max(maxAbs * 1.05f, 1e-4f);

        _vertexColorMat = Resources.Load<Material>("Materials/VertexColorUnlit");
        if (_vertexColorMat == null)
        {
            Shader s = Shader.Find("Custom/VertexColorUnlit") ?? Shader.Find("Universal Render Pipeline/Unlit");
            _vertexColorMat = new Material(s);
        }

        _dataLoaded = true;
    }

    private Vector3 MapPoint(int e, int t)
    {
        int idx = e * _t + t;
        float x = (_wr[idx] / _rmax) * (BoxW * 0.5f);
        float z = (_wi[idx] / _rmax) * (BoxW * 0.5f);
        float logV = Mathf.Log10(Mathf.Max(_fz[t], _floor));
        float yNorm = Mathf.InverseLerp(MinLog, MaxLog, logV);
        float y = (yNorm - 0.5f) * BoxH;
        return new Vector3(x, y, z);
    }

    private void BuildScene3D()
    {
        _sceneRoot = new GameObject("Inflation3DScene");

        BuildGridBox();
        BuildFiberMesh();
        BuildMarkerMesh();
        BuildAxisLabels();
    }

    private void BuildGridBox()
    {
        var go = new GameObject("GridBox");
        go.transform.SetParent(_sceneRoot.transform, false);
        var mf = go.AddComponent<MeshFilter>();
        var mr = go.AddComponent<MeshRenderer>();
        mr.sharedMaterial = _vertexColorMat;

        float hw = BoxW * 0.5f;
        float hh = BoxH * 0.5f;

        var verts = new System.Collections.Generic.List<Vector3>();
        var cols = new System.Collections.Generic.List<Color>();
        var indices = new System.Collections.Generic.List<int>();

        void AddSeg(Vector3 a, Vector3 b, Color c)
        {
            int idx = verts.Count;
            verts.Add(a); verts.Add(b);
            cols.Add(c); cols.Add(c);
            indices.Add(idx); indices.Add(idx + 1);
        }

        Color edgeCol = new Color(0.25f, 0.30f, 0.38f, 0.95f);
        Color gridCol = new Color(0.65f, 0.72f, 0.80f, 0.65f);

        // Outer box 12 edges
        Vector3[] c =
        {
            new Vector3(-hw, -hh, -hw), new Vector3( hw, -hh, -hw),
            new Vector3( hw, -hh,  hw), new Vector3(-hw, -hh,  hw),
            new Vector3(-hw,  hh, -hw), new Vector3( hw,  hh, -hw),
            new Vector3( hw,  hh,  hw), new Vector3(-hw,  hh,  hw),
        };
        for (int i = 0; i < 4; i++)
        {
            AddSeg(c[i], c[(i + 1) % 4], edgeCol);
            AddSeg(c[i + 4], c[((i + 1) % 4) + 4], edgeCol);
            AddSeg(c[i], c[i + 4], edgeCol);
        }

        // Bottom & top floor grid (5x5)
        for (int i = 1; i < 6; i++)
        {
            float d = -hw + (BoxW * i) / 6f;
            AddSeg(new Vector3(d, -hh, -hw), new Vector3(d, -hh, hw), gridCol);
            AddSeg(new Vector3(-hw, -hh, d), new Vector3(hw, -hh, d), gridCol);
        }

        // Horizontal Y-level grid rings at log10 = -5, -10, -15, -20, -25, -30
        float[] yTicks = { -5f, -10f, -15f, -20f, -25f, -30f };
        foreach (float lv in yTicks)
        {
            float y = (Mathf.InverseLerp(MinLog, MaxLog, lv) - 0.5f) * BoxH;
            AddSeg(new Vector3(-hw, y, -hw), new Vector3( hw, y, -hw), gridCol);
            AddSeg(new Vector3( hw, y, -hw), new Vector3( hw, y,  hw), gridCol);
            AddSeg(new Vector3( hw, y,  hw), new Vector3(-hw, y,  hw), gridCol);
            AddSeg(new Vector3(-hw, y,  hw), new Vector3(-hw, y, -hw), gridCol);
        }

        var mesh = new Mesh { name = "GridBoxMesh" };
        mesh.SetVertices(verts);
        mesh.SetColors(cols);
        mesh.SetIndices(indices, MeshTopology.Lines, 0);
        mesh.RecalculateBounds();
        mf.sharedMesh = mesh;
    }

    private void BuildFiberMesh()
    {
        var go = new GameObject("Fibers");
        go.transform.SetParent(_sceneRoot.transform, false);
        var mf = go.AddComponent<MeshFilter>();
        var mr = go.AddComponent<MeshRenderer>();
        mr.sharedMaterial = _vertexColorMat;

        int vertCount = _t * _m;
        var verts = new Vector3[vertCount];
        var cols = new Color[vertCount];

        for (int t = 0; t < _t; t++)
        {
            for (int e = 0; e < _m; e++)
            {
                int vi = t * _m + e;
                verts[vi] = MapPoint(e, t);
                Color col = Palette[e % Palette.Length];
                col.a = 0.90f;
                cols[vi] = col;
            }
        }

        int segCount = (_t - 1) * _m;
        var indices = new int[segCount * 2];
        int ptr = 0;
        for (int t = 0; t < _t - 1; t++)
        {
            int base0 = t * _m;
            int base1 = (t + 1) * _m;
            for (int e = 0; e < _m; e++)
            {
                indices[ptr++] = base0 + e;
                indices[ptr++] = base1 + e;
            }
        }

        _fiberMesh = new Mesh
        {
            name = "FibersMesh",
            indexFormat = IndexFormat.UInt32
        };
        _fiberMesh.SetVertices(verts);
        _fiberMesh.SetColors(cols);
        _fiberMesh.SetIndices(indices, MeshTopology.Lines, 0);
        _fiberMesh.bounds = new Bounds(Vector3.zero, new Vector3(BoxW, BoxH, BoxW));
        mf.sharedMesh = _fiberMesh;
    }

    private void BuildMarkerMesh()
    {
        var go = new GameObject("CurrentMarkers");
        go.transform.SetParent(_sceneRoot.transform, false);
        var mf = go.AddComponent<MeshFilter>();
        var mr = go.AddComponent<MeshRenderer>();
        mr.sharedMaterial = _vertexColorMat;

        _markerVerts = new Vector3[_m * 4];
        var cols = new Color[_m * 4];
        var tris = new int[_m * 6];

        for (int e = 0; e < _m; e++)
        {
            Color col = Palette[e % Palette.Length];
            col.a = 1f;
            int vBase = e * 4;
            cols[vBase] = col;
            cols[vBase + 1] = col;
            cols[vBase + 2] = col;
            cols[vBase + 3] = col;

            int tBase = e * 6;
            tris[tBase] = vBase;
            tris[tBase + 1] = vBase + 1;
            tris[tBase + 2] = vBase + 2;
            tris[tBase + 3] = vBase;
            tris[tBase + 4] = vBase + 2;
            tris[tBase + 5] = vBase + 3;
        }

        _markerMesh = new Mesh { name = "MarkersMesh" };
        _markerMesh.vertices = _markerVerts;
        _markerMesh.colors = cols;
        _markerMesh.triangles = tris;
        _markerMesh.bounds = new Bounds(Vector3.zero, new Vector3(BoxW, BoxH, BoxW));
        mf.sharedMesh = _markerMesh;
    }

    private void UpdateMarkersBillboard()
    {
        if (_markerMesh == null || _camera == null) return;
        const float halfSize = 0.065f;
        Vector3 right = _camera.transform.right * halfSize;
        Vector3 up = _camera.transform.up * halfSize;

        for (int e = 0; e < _m; e++)
        {
            Vector3 p = MapPoint(e, _currentStep);
            int vBase = e * 4;
            _markerVerts[vBase] = p - right - up;
            _markerVerts[vBase + 1] = p - right + up;
            _markerVerts[vBase + 2] = p + right + up;
            _markerVerts[vBase + 3] = p + right - up;
        }
        _markerMesh.vertices = _markerVerts;
    }

    private void BuildAxisLabels()
    {
        float hw = BoxW * 0.5f;
        float hh = BoxH * 0.5f;
        _billboardCount = 0;

        AddWorldLabel("Re w", new Vector3(0f, -hh - 0.45f, -hw - 0.35f), 2.8f);
        AddWorldLabel("Im w", new Vector3(-hw - 0.45f, -hh - 0.45f, 0f), 2.8f);
        AddWorldLabel("log₁₀(H⊥/H) ← 膨張", new Vector3(-hw - 1.15f, 0f, hw + 0.35f), 2.4f);

        float[] yTicks = { 0f, -5f, -10f, -15f, -20f, -25f, -30f };
        foreach (float lv in yTicks)
        {
            float y = (Mathf.InverseLerp(MinLog, MaxLog, lv) - 0.5f) * BoxH;
            AddWorldLabel(((int)lv).ToString(), new Vector3(-hw - 0.38f, y, hw + 0.15f), 2.3f);
        }
    }

    private void AddWorldLabel(string text, Vector3 pos, float fontSize)
    {
        var go = new GameObject("AxisLabel_" + text);
        go.transform.SetParent(_sceneRoot.transform, false);
        go.transform.localPosition = pos;
        var tmp = go.AddComponent<TextMeshPro>();
        if (_font != null) tmp.font = _font;
        tmp.text = text;
        tmp.fontSize = fontSize;
        tmp.color = new Color(0.15f, 0.2f, 0.28f, 1f);
        tmp.alignment = TextAlignmentOptions.Center;
        tmp.rectTransform.sizeDelta = new Vector2(8f, 1.2f);
        if (_billboardCount < _billboards.Length)
        {
            _billboards[_billboardCount++] = go.transform;
        }
    }

    private void BuildUI(RectTransform panelRoot)
    {
        var containerGo = new GameObject("InflationUI", typeof(RectTransform));
        _uiContainer = containerGo.GetComponent<RectTransform>();
        _uiContainer.SetParent(panelRoot, false);
        _uiContainer.anchorMin = Vector2.zero;
        _uiContainer.anchorMax = Vector2.one;
        _uiContainer.offsetMin = Vector2.zero;
        _uiContainer.offsetMax = Vector2.zero;

        // Top title overlay inside panel
        var titleText = CreateTMP(_uiContainer, "Title",
            "関係波の膨張＝インフレーション（N=40・M=780本の波）", 20, TextAlignmentOptions.Center,
            new Color(0.10f, 0.14f, 0.22f, 1f));
        var titleRt = titleText.rectTransform;
        titleRt.anchorMin = new Vector2(0f, 1f);
        titleRt.anchorMax = new Vector2(1f, 1f);
        titleRt.pivot = new Vector2(0.5f, 1f);
        titleRt.sizeDelta = new Vector2(0f, 36f);
        titleRt.anchoredPosition = new Vector2(0f, -6f);

        // Bottom control + 2D curve dock
        var dockGo = new GameObject("BottomDock", typeof(RectTransform), typeof(Image));
        var dockRt = dockGo.GetComponent<RectTransform>();
        dockRt.SetParent(_uiContainer, false);
        dockRt.anchorMin = new Vector2(0f, 0f);
        dockRt.anchorMax = new Vector2(1f, 0f);
        dockRt.pivot = new Vector2(0.5f, 0f);
        dockRt.sizeDelta = new Vector2(0f, 154f);
        dockRt.anchoredPosition = Vector2.zero;
        dockGo.GetComponent<Image>().color = new Color(1f, 1f, 1f, 0.92f);

        // Controls row (top 44px of dock)
        var ctlGo = new GameObject("ControlsRow", typeof(RectTransform));
        var ctlRt = ctlGo.GetComponent<RectTransform>();
        ctlRt.SetParent(dockRt, false);
        ctlRt.anchorMin = new Vector2(0f, 1f);
        ctlRt.anchorMax = new Vector2(1f, 1f);
        ctlRt.pivot = new Vector2(0.5f, 1f);
        ctlRt.sizeDelta = new Vector2(-24f, 40f);
        ctlRt.anchoredPosition = new Vector2(0f, -6f);

        // Play / Pause Button
        var btnGo = new GameObject("PlayButton", typeof(RectTransform), typeof(Image), typeof(Button));
        var btnRt = btnGo.GetComponent<RectTransform>();
        btnRt.SetParent(ctlRt, false);
        btnRt.anchorMin = new Vector2(0f, 0.5f);
        btnRt.anchorMax = new Vector2(0f, 0.5f);
        btnRt.pivot = new Vector2(0f, 0.5f);
        btnRt.sizeDelta = new Vector2(96f, 34f);
        btnRt.anchoredPosition = Vector2.zero;
        btnGo.GetComponent<Image>().color = new Color(0.18f, 0.45f, 0.92f, 1f);
        var btn = btnGo.GetComponent<Button>();
        btn.onClick.AddListener(TogglePlay);

        _playBtnLabel = CreateTMP(btnRt, "Label", "▶ 再生", 15, TextAlignmentOptions.Center, Color.white);
        StretchFull(_playBtnLabel.rectTransform);

        // Status label (right side)
        _statusLabel = CreateTMP(ctlRt, "Status", "step 0/500  H⊥/H = 1.00e-31", 15,
            TextAlignmentOptions.Right, new Color(0.15f, 0.2f, 0.28f, 1f));
        var stRt = _statusLabel.rectTransform;
        stRt.anchorMin = new Vector2(1f, 0.5f);
        stRt.anchorMax = new Vector2(1f, 0.5f);
        stRt.pivot = new Vector2(1f, 0.5f);
        stRt.sizeDelta = new Vector2(250f, 34f);
        stRt.anchoredPosition = Vector2.zero;

        // Slider (between button and status)
        _slider = CreateSlider(ctlRt, 108f, 262f);
        _slider.wholeNumbers = true;
        _slider.minValue = 0;
        _slider.maxValue = _t - 1;
        _slider.onValueChanged.AddListener(val =>
        {
            if (_suppressSliderCallback) return;
            SetPlaying(false);
            SetStep(Mathf.RoundToInt(val));
        });

        // 2D Inflation Curve container
        var graphBgGo = new GameObject("CurveArea", typeof(RectTransform), typeof(Image));
        var graphBgRt = graphBgGo.GetComponent<RectTransform>();
        graphBgRt.SetParent(dockRt, false);
        graphBgRt.anchorMin = new Vector2(0f, 0f);
        graphBgRt.anchorMax = new Vector2(1f, 1f);
        graphBgRt.offsetMin = new Vector2(64f, 10f);
        graphBgRt.offsetMax = new Vector2(-16f, -50f);
        graphBgGo.GetComponent<Image>().color = new Color(0.95f, 0.97f, 0.99f, 1f);

        var curveGo = new GameObject("InflationCurve", typeof(RectTransform), typeof(InflationCurveGraphic));
        var curveRt = curveGo.GetComponent<RectTransform>();
        curveRt.SetParent(graphBgRt, false);
        StretchFull(curveRt);
        curveRt.offsetMin = new Vector2(4f, 4f);
        curveRt.offsetMax = new Vector2(-4f, -4f);
        _curveGraphic = curveGo.GetComponent<InflationCurveGraphic>();
        _curveGraphic.SetData(_fz, _floor);

        // Y-axis caption on left of graph
        var yCap = CreateTMP(dockRt, "CurveYLabel", "H⊥/H\n10^0\n\n10^-31", 11,
            TextAlignmentOptions.Right, new Color(0.35f, 0.4f, 0.48f, 1f));
        var yCapRt = yCap.rectTransform;
        yCapRt.anchorMin = new Vector2(0f, 0f);
        yCapRt.anchorMax = new Vector2(0f, 1f);
        yCapRt.offsetMin = new Vector2(4f, 8f);
        yCapRt.offsetMax = new Vector2(58f, -50f);
    }

    private void TogglePlay()
    {
        if (!_isPlaying && _currentStep >= _t - 1)
        {
            SetStep(0);
        }
        SetPlaying(!_isPlaying);
    }

    private void SetPlaying(bool play)
    {
        _isPlaying = play;
        _playTimer = 0f;
        if (_playBtnLabel != null)
        {
            _playBtnLabel.text = play ? "⏸ 停止" : "▶ 再生";
        }
    }

    private void SetStep(int step)
    {
        _currentStep = Mathf.Clamp(step, 0, _t - 1);

        if (_fiberMesh != null)
        {
            int indexCount = _currentStep * _m * 2;
            _fiberMesh.SetSubMesh(0, new SubMeshDescriptor(0, indexCount, MeshTopology.Lines));
        }

        UpdateMarkersBillboard();

        if (_curveGraphic != null) _curveGraphic.SetStep(_currentStep);

        if (_slider != null)
        {
            _suppressSliderCallback = true;
            _slider.value = _currentStep;
            _suppressSliderCallback = false;
        }

        if (_statusLabel != null)
        {
            _statusLabel.text = $"step {_currentStep}/{_t - 1}  H⊥/H = {_fz[_currentStep]:0.00e+00}";
        }
    }

    private TextMeshProUGUI CreateTMP(Transform parent, string name, string text, float size,
        TextAlignmentOptions align, Color color)
    {
        var go = new GameObject(name, typeof(RectTransform), typeof(TextMeshProUGUI));
        go.transform.SetParent(parent, false);
        var tmp = go.GetComponent<TextMeshProUGUI>();
        if (_font != null) tmp.font = _font;
        tmp.text = text;
        tmp.fontSize = size;
        tmp.alignment = align;
        tmp.color = color;
        return tmp;
    }

    private static void StretchFull(RectTransform rt)
    {
        rt.anchorMin = Vector2.zero;
        rt.anchorMax = Vector2.one;
        rt.offsetMin = Vector2.zero;
        rt.offsetMax = Vector2.zero;
    }

    private static Slider CreateSlider(Transform parent, float leftPad, float rightPad)
    {
        var rootGo = new GameObject("StepSlider", typeof(RectTransform), typeof(Slider));
        var rootRt = rootGo.GetComponent<RectTransform>();
        rootRt.SetParent(parent, false);
        rootRt.anchorMin = new Vector2(0f, 0.5f);
        rootRt.anchorMax = new Vector2(1f, 0.5f);
        rootRt.offsetMin = new Vector2(leftPad, -14f);
        rootRt.offsetMax = new Vector2(-rightPad, 14f);

        var bgGo = new GameObject("Background", typeof(RectTransform), typeof(Image));
        var bgRt = bgGo.GetComponent<RectTransform>();
        bgRt.SetParent(rootRt, false);
        bgRt.anchorMin = new Vector2(0f, 0.35f);
        bgRt.anchorMax = new Vector2(1f, 0.65f);
        bgRt.offsetMin = Vector2.zero;
        bgRt.offsetMax = Vector2.zero;
        bgGo.GetComponent<Image>().color = new Color(0.82f, 0.85f, 0.90f, 1f);

        var fillArea = new GameObject("Fill Area", typeof(RectTransform));
        var faRt = fillArea.GetComponent<RectTransform>();
        faRt.SetParent(rootRt, false);
        faRt.anchorMin = new Vector2(0f, 0.35f);
        faRt.anchorMax = new Vector2(1f, 0.65f);
        faRt.offsetMin = Vector2.zero;
        faRt.offsetMax = Vector2.zero;

        var fillGo = new GameObject("Fill", typeof(RectTransform), typeof(Image));
        var fillRt = fillGo.GetComponent<RectTransform>();
        fillRt.SetParent(faRt, false);
        StretchFull(fillRt);
        fillGo.GetComponent<Image>().color = new Color(0.23f, 0.51f, 0.96f, 1f);

        var handleArea = new GameObject("Handle Slide Area", typeof(RectTransform));
        var haRt = handleArea.GetComponent<RectTransform>();
        haRt.SetParent(rootRt, false);
        StretchFull(haRt);

        var handleGo = new GameObject("Handle", typeof(RectTransform), typeof(Image));
        var handleRt = handleGo.GetComponent<RectTransform>();
        handleRt.SetParent(haRt, false);
        handleRt.sizeDelta = new Vector2(20f, 20f);
        handleGo.GetComponent<Image>().color = new Color(0.12f, 0.36f, 0.82f, 1f);

        var slider = rootGo.GetComponent<Slider>();
        slider.fillRect = fillRt;
        slider.handleRect = handleRt;
        slider.targetGraphic = handleGo.GetComponent<Image>();
        slider.direction = Slider.Direction.LeftToRight;
        return slider;
    }

    private static Color HexColor(int rgb)
    {
        return new Color(
            ((rgb >> 16) & 0xFF) / 255f,
            ((rgb >> 8) & 0xFF) / 255f,
            (rgb & 0xFF) / 255f,
            1f).linear;
    }
}
