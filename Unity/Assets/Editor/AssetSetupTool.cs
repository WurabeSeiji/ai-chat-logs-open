using System.IO;
using TMPro;
using UnityEditor;
using UnityEngine;
using UnityEngine.TextCore.LowLevel;

public static class AssetSetupTool
{
    private const string MaterialPath = "Assets/Resources/Materials/VertexColorUnlit.mat";
    private const string FontAssetPath = "Assets/Resources/Fonts/JapaneseSDF.asset";
    private const string SourceFontPath = "/System/Library/Fonts/Supplemental/Arial Unicode.ttf";

    [MenuItem("Tools/Setup Runtime Assets")]
    public static void SetupRuntimeAssets()
    {
        EnsureFolder("Assets/Resources/Materials");
        EnsureFolder("Assets/Resources/Fonts");

        EnsureTMPEssentialResources();
        SetupVertexColorMaterial();
        SetupJapaneseFontAsset();

        AssetDatabase.SaveAssets();
        AssetDatabase.Refresh();
        Debug.Log("[AssetSetupTool] Runtime assets setup completed.");

        if (Application.isBatchMode)
        {
            EditorApplication.Exit(0);
        }
    }

    private static void EnsureTMPEssentialResources()
    {
        if (Resources.Load<TMP_Settings>("TMP Settings") != null)
        {
            return;
        }
        string pkgPath = Path.GetFullPath("Packages/com.unity.ugui/Package Resources/TMP Essential Resources.unitypackage");
        if (File.Exists(pkgPath))
        {
            AssetDatabase.ImportPackage(pkgPath, false);
            AssetDatabase.Refresh();
            Debug.Log("[AssetSetupTool] Imported TMP Essential Resources.");
        }
        else
        {
            Debug.LogWarning($"[AssetSetupTool] TMP Essential Resources package not found: {pkgPath}");
        }
    }

    private static void SetupVertexColorMaterial()
    {
        Shader shader = Shader.Find("Custom/VertexColorUnlit");
        if (shader == null)
        {
            Debug.LogError("[AssetSetupTool] Shader 'Custom/VertexColorUnlit' not found.");
            return;
        }

        Material mat = AssetDatabase.LoadAssetAtPath<Material>(MaterialPath);
        if (mat == null)
        {
            mat = new Material(shader);
            AssetDatabase.CreateAsset(mat, MaterialPath);
        }
        else
        {
            mat.shader = shader;
            EditorUtility.SetDirty(mat);
        }
    }

    private static void SetupJapaneseFontAsset()
    {
        if (!File.Exists(SourceFontPath))
        {
            Debug.LogWarning($"[AssetSetupTool] Source font not found: {SourceFontPath}");
            return;
        }

        if (AssetDatabase.LoadAssetAtPath<TMP_FontAsset>(FontAssetPath) != null)
        {
            AssetDatabase.DeleteAsset(FontAssetPath);
        }

        Font sourceFont = new Font(SourceFontPath);
        TMP_FontAsset fontAsset = TMP_FontAsset.CreateFontAsset(
            sourceFont,
            36,
            4,
            GlyphRenderMode.SDFAA,
            1024,
            1024,
            AtlasPopulationMode.Dynamic,
            false
        );

        if (fontAsset == null)
        {
            Debug.LogError("[AssetSetupTool] Failed to create TMP_FontAsset.");
            return;
        }

        fontAsset.name = "JapaneseSDF";
        string charset = BuildCharset();
        fontAsset.TryAddCharacters(charset, out string missing);
        if (!string.IsNullOrEmpty(missing))
        {
            Debug.LogWarning($"[AssetSetupTool] Missing characters: {missing}");
        }

        fontAsset.atlasPopulationMode = AtlasPopulationMode.Static;
        fontAsset.ReadFontAssetDefinition();

        AssetDatabase.CreateAsset(fontAsset, FontAssetPath);
        if (fontAsset.atlasTexture != null)
        {
            fontAsset.atlasTexture.name = "JapaneseSDF Atlas";
            AssetDatabase.AddObjectToAsset(fontAsset.atlasTexture, fontAsset);
        }
        if (fontAsset.material != null)
        {
            fontAsset.material.name = "JapaneseSDF Material";
            AssetDatabase.AddObjectToAsset(fontAsset.material, fontAsset);
        }
        EditorUtility.SetDirty(fontAsset);
    }

    private static string BuildCharset()
    {
        var sb = new System.Text.StringBuilder();
        // Basic ASCII
        for (int c = 32; c <= 126; c++) sb.Append((char)c);
        // Full-width / Japanese punctuation & symbols
        sb.Append("、。・：；！？（）「」『』【】＝＋－×÷←→↑↓⊥₁₀・〜…▶⏸■□●○");
        // Hiragana
        for (int c = 0x3041; c <= 0x3096; c++) sb.Append((char)c);
        // Katakana
        for (int c = 0x30A1; c <= 0x30FC; c++) sb.Append((char)c);
        // Kanji used in physics / UI / tabs
        sb.Append("関係波膨張本位相複素平面共回転空間軸縦横時間再生一時停止停止速度初期戻一覧選択画面切替相対独立曲線追試実験次元生成構造電子反跳光子");
        return sb.ToString();
    }

    private static void EnsureFolder(string path)
    {
        if (AssetDatabase.IsValidFolder(path)) return;
        string parent = Path.GetDirectoryName(path)?.Replace("\\", "/");
        string name = Path.GetFileName(path);
        if (!string.IsNullOrEmpty(parent) && !AssetDatabase.IsValidFolder(parent))
        {
            EnsureFolder(parent);
        }
        AssetDatabase.CreateFolder(parent, name);
    }
}
