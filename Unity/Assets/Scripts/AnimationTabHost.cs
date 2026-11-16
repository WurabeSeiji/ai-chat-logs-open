using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.InputSystem.UI;
using UnityEngine.UI;

public class AnimationTabHost : MonoBehaviour
{
    private const float TabBarHeight = 52f;

    private static readonly Color TabBarBg = new Color(0.11f, 0.14f, 0.20f, 1f);
    private static readonly Color TabActiveBg = new Color(0.23f, 0.51f, 0.96f, 1f);
    private static readonly Color TabInactiveBg = new Color(0.19f, 0.23f, 0.31f, 1f);

    private readonly List<IAnimationModule> _modules = new List<IAnimationModule>();
    private readonly List<Image> _tabImages = new List<Image>();
    private readonly List<TextMeshProUGUI> _tabLabels = new List<TextMeshProUGUI>();

    private RectTransform _tabContent;
    private RectTransform _panelRoot;
    private Camera _mainCamera;
    private TMP_FontAsset _font;
    private int _activeIndex = -1;

    [RuntimeInitializeOnLoadMethod(RuntimeInitializeLoadType.AfterSceneLoad)]
    private static void Bootstrap()
    {
        if (FindFirstObjectByType<AnimationTabHost>() != null) return;
        var hostGo = new GameObject("AnimationTabHost");
        DontDestroyOnLoad(hostGo);
        hostGo.AddComponent<AnimationTabHost>();
    }

    private void Start()
    {
        _font = Resources.Load<TMP_FontAsset>("Fonts/JapaneseSDF");
        _mainCamera = Camera.main;
        if (_mainCamera == null)
        {
            _mainCamera = FindFirstObjectByType<Camera>();
        }

        EnsureEventSystem();
        BuildShellUI();

        // Register animation modules here (add future animations to this list)
        RegisterModule(gameObject.AddComponent<InflationAnimationModule>());

        if (_modules.Count > 0)
        {
            SelectTab(0);
        }
    }

    public void RegisterModule(IAnimationModule module)
    {
        int index = _modules.Count;
        _modules.Add(module);
        CreateTabButton(index, module.Title);
    }

    public void SelectTab(int index)
    {
        if (index < 0 || index >= _modules.Count) return;

        if (_activeIndex >= 0 && _activeIndex < _modules.Count)
        {
            _modules[_activeIndex].Deactivate();
        }

        _activeIndex = index;
        for (int i = 0; i < _tabImages.Count; i++)
        {
            bool active = (i == index);
            _tabImages[i].color = active ? TabActiveBg : TabInactiveBg;
            _tabLabels[i].color = active ? Color.white : new Color(0.80f, 0.84f, 0.90f, 1f);
        }

        _modules[index].Activate(_panelRoot, _mainCamera, _font);
    }

    private void EnsureEventSystem()
    {
        if (FindFirstObjectByType<EventSystem>() != null) return;
        var esGo = new GameObject("EventSystem", typeof(EventSystem), typeof(InputSystemUIInputModule));
        DontDestroyOnLoad(esGo);
    }

    private void BuildShellUI()
    {
        var canvasGo = new GameObject("MainUICanvas", typeof(RectTransform), typeof(Canvas), typeof(CanvasScaler), typeof(GraphicRaycaster));
        canvasGo.transform.SetParent(transform, false);

        var canvas = canvasGo.GetComponent<Canvas>();
        canvas.renderMode = RenderMode.ScreenSpaceOverlay;
        canvas.sortingOrder = 10;

        var scaler = canvasGo.GetComponent<CanvasScaler>();
        scaler.uiScaleMode = CanvasScaler.ScaleMode.ScaleWithScreenSize;
        scaler.referenceResolution = new Vector2(1280f, 720f);
        scaler.matchWidthOrHeight = 0.5f;

        // Top tab bar with horizontal ScrollRect
        var tabBarGo = new GameObject("TopTabBar", typeof(RectTransform), typeof(Image), typeof(ScrollRect));
        var tabBarRt = tabBarGo.GetComponent<RectTransform>();
        tabBarRt.SetParent(canvasGo.transform, false);
        tabBarRt.anchorMin = new Vector2(0f, 1f);
        tabBarRt.anchorMax = new Vector2(1f, 1f);
        tabBarRt.pivot = new Vector2(0.5f, 1f);
        tabBarRt.sizeDelta = new Vector2(0f, TabBarHeight);
        tabBarRt.anchoredPosition = Vector2.zero;
        tabBarGo.GetComponent<Image>().color = TabBarBg;

        var viewportGo = new GameObject("Viewport", typeof(RectTransform), typeof(RectMask2D));
        var vpRt = viewportGo.GetComponent<RectTransform>();
        vpRt.SetParent(tabBarRt, false);
        vpRt.anchorMin = Vector2.zero;
        vpRt.anchorMax = Vector2.one;
        vpRt.offsetMin = new Vector2(8f, 4f);
        vpRt.offsetMax = new Vector2(-8f, -4f);

        var contentGo = new GameObject("TabContent", typeof(RectTransform), typeof(HorizontalLayoutGroup), typeof(ContentSizeFitter));
        _tabContent = contentGo.GetComponent<RectTransform>();
        _tabContent.SetParent(vpRt, false);
        _tabContent.anchorMin = new Vector2(0f, 0f);
        _tabContent.anchorMax = new Vector2(0f, 1f);
        _tabContent.pivot = new Vector2(0f, 0.5f);
        _tabContent.anchoredPosition = Vector2.zero;

        var hlg = contentGo.GetComponent<HorizontalLayoutGroup>();
        hlg.spacing = 8f;
        hlg.childAlignment = TextAnchor.MiddleLeft;
        hlg.childControlWidth = false;
        hlg.childControlHeight = true;
        hlg.childForceExpandWidth = false;
        hlg.childForceExpandHeight = true;

        var csf = contentGo.GetComponent<ContentSizeFitter>();
        csf.horizontalFit = ContentSizeFitter.FitMode.PreferredSize;
        csf.verticalFit = ContentSizeFitter.FitMode.Unconstrained;

        var scroll = tabBarGo.GetComponent<ScrollRect>();
        scroll.viewport = vpRt;
        scroll.content = _tabContent;
        scroll.horizontal = true;
        scroll.vertical = false;
        scroll.movementType = ScrollRect.MovementType.Clamped;

        // Lower content panel (everything below TopTabBar)
        var panelGo = new GameObject("AnimationPanel", typeof(RectTransform));
        _panelRoot = panelGo.GetComponent<RectTransform>();
        _panelRoot.SetParent(canvasGo.transform, false);
        _panelRoot.anchorMin = Vector2.zero;
        _panelRoot.anchorMax = Vector2.one;
        _panelRoot.offsetMin = Vector2.zero;
        _panelRoot.offsetMax = new Vector2(0f, -TabBarHeight);
    }

    private void CreateTabButton(int index, string title)
    {
        var btnGo = new GameObject("Tab_" + title, typeof(RectTransform), typeof(Image), typeof(Button), typeof(LayoutElement));
        var rt = btnGo.GetComponent<RectTransform>();
        rt.SetParent(_tabContent, false);

        float width = Mathf.Max(176f, title.Length * 22f + 48f);
        rt.sizeDelta = new Vector2(width, TabBarHeight - 8f);
        var le = btnGo.GetComponent<LayoutElement>();
        le.preferredWidth = width;
        le.minWidth = width;

        var img = btnGo.GetComponent<Image>();
        img.color = TabInactiveBg;
        _tabImages.Add(img);

        var lblGo = new GameObject("Label", typeof(RectTransform), typeof(TextMeshProUGUI));
        var lblRt = lblGo.GetComponent<RectTransform>();
        lblRt.SetParent(rt, false);
        lblRt.anchorMin = Vector2.zero;
        lblRt.anchorMax = Vector2.one;
        lblRt.offsetMin = new Vector2(10f, 2f);
        lblRt.offsetMax = new Vector2(-10f, -2f);

        var tmp = lblGo.GetComponent<TextMeshProUGUI>();
        if (_font != null) tmp.font = _font;
        tmp.text = title;
        tmp.fontSize = 18f;
        tmp.textWrappingMode = TextWrappingModes.NoWrap;
        tmp.overflowMode = TextOverflowModes.Overflow;
        tmp.alignment = TextAlignmentOptions.Center;
        tmp.color = Color.white;
        _tabLabels.Add(tmp);

        int captured = index;
        btnGo.GetComponent<Button>().onClick.AddListener(() => SelectTab(captured));
    }
}
