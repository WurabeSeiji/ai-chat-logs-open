using UnityEngine;
using UnityEngine.UI;

[RequireComponent(typeof(CanvasRenderer))]
public class InflationCurveGraphic : MaskableGraphic
{
    private float[] _fz;
    private float _minLog = -31f;
    private float _maxLog = 0f;
    private int _currentStep;

    private static readonly Color CurveColor = new Color(0.173f, 0.627f, 0.173f, 1f); // #2ca02c
    private static readonly Color MarkerColor = new Color(1.0f, 0.498f, 0.055f, 1f); // #ff7f0e
    private static readonly Color GuideColor = new Color(0.53f, 0.53f, 0.53f, 0.85f);
    private static readonly Color GridColor = new Color(0.85f, 0.88f, 0.92f, 1f);

    public void SetData(float[] fz, float floor)
    {
        _fz = fz;
        _minLog = Mathf.Log10(Mathf.Max(floor, 1e-35f));
        _maxLog = 0f;
        SetVerticesDirty();
    }

    public void SetStep(int step)
    {
        _currentStep = step;
        SetVerticesDirty();
    }

    protected override void OnPopulateMesh(VertexHelper vh)
    {
        vh.Clear();
        Rect r = GetPixelAdjustedRect();
        if (r.width <= 2f || r.height <= 2f) return;

        // Subtle horizontal grid lines at log10 = 0, -10, -20, -30
        float[] ticks = { 0f, -10f, -20f, -30f };
        foreach (float t in ticks)
        {
            float y = MapY(t, r);
            AddLine(vh, new Vector2(r.xMin, y), new Vector2(r.xMax, y), 1f, GridColor);
        }

        if (_fz == null || _fz.Length < 2) return;

        int n = _fz.Length;
        Vector2 prev = Vector2.zero;
        for (int i = 0; i < n; i++)
        {
            float x = r.xMin + (r.width * i) / (n - 1);
            float logV = Mathf.Log10(Mathf.Max(_fz[i], 1e-35f));
            float y = MapY(logV, r);
            Vector2 curr = new Vector2(x, y);
            if (i > 0)
            {
                AddLine(vh, prev, curr, 2.2f, CurveColor);
            }
            prev = curr;
        }

        int k = Mathf.Clamp(_currentStep, 0, n - 1);
        float kx = r.xMin + (r.width * k) / (n - 1);
        float kLog = Mathf.Log10(Mathf.Max(_fz[k], 1e-35f));
        float ky = MapY(kLog, r);

        // Dashed vertical guide from floor to current point
        const float dashLen = 4f;
        const float gapLen = 3f;
        for (float y = r.yMin; y < ky; y += dashLen + gapLen)
        {
            float yEnd = Mathf.Min(y + dashLen, ky);
            AddLine(vh, new Vector2(kx, y), new Vector2(kx, yEnd), 1.5f, GuideColor);
        }

        // Current step orange dot
        AddCircle(vh, new Vector2(kx, ky), 5.0f, 12, MarkerColor);
    }

    private float MapY(float logVal, Rect r)
    {
        float t = Mathf.InverseLerp(_minLog, _maxLog, logVal);
        return r.yMin + t * r.height;
    }

    private static void AddLine(VertexHelper vh, Vector2 a, Vector2 b, float thickness, Color col)
    {
        Vector2 dir = b - a;
        float len = dir.magnitude;
        if (len < 1e-4f) return;
        Vector2 normal = new Vector2(-dir.y, dir.x) * (thickness * 0.5f / len);

        int idx = vh.currentVertCount;
        vh.AddVert(a - normal, col, Vector4.zero);
        vh.AddVert(a + normal, col, Vector4.zero);
        vh.AddVert(b + normal, col, Vector4.zero);
        vh.AddVert(b - normal, col, Vector4.zero);

        vh.AddTriangle(idx, idx + 1, idx + 2);
        vh.AddTriangle(idx, idx + 2, idx + 3);
    }

    private static void AddCircle(VertexHelper vh, Vector2 center, float radius, int segments, Color col)
    {
        int centerIdx = vh.currentVertCount;
        vh.AddVert(center, col, Vector4.zero);
        for (int i = 0; i <= segments; i++)
        {
            float ang = (Mathf.PI * 2f * i) / segments;
            Vector2 p = center + new Vector2(Mathf.Cos(ang), Mathf.Sin(ang)) * radius;
            vh.AddVert(p, col, Vector4.zero);
        }
        for (int i = 0; i < segments; i++)
        {
            vh.AddTriangle(centerIdx, centerIdx + 1 + i, centerIdx + 2 + i);
        }
    }
}
