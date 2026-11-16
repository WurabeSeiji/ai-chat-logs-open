using System;
using System.IO;
using System.Linq;
using UnityEditor;
using UnityEditor.Build.Reporting;
using UnityEngine;

public static class WebGLBuilder
{
    private const string OutputPath = "Builds/WebGL";

    [MenuItem("Build/Build WebGL (WASM)")]
    public static void BuildWebGL()
    {
        // PC・スマホの静的ホスティング（Apache / GitHub Pages / python http.server 等）で確実に起動する設定
        PlayerSettings.WebGL.compressionFormat = WebGLCompressionFormat.Gzip;
        PlayerSettings.WebGL.decompressionFallback = true;
        PlayerSettings.WebGL.threadsSupport = false;
        PlayerSettings.runInBackground = true;

        string[] scenes = EditorBuildSettings.scenes
            .Where(s => s.enabled)
            .Select(s => s.path)
            .ToArray();

        if (scenes.Length == 0)
        {
            scenes = new[] { "Assets/Scenes/SampleScene.unity" };
        }

        var options = new BuildPlayerOptions
        {
            scenes = scenes,
            locationPathName = OutputPath,
            target = BuildTarget.WebGL,
            options = BuildOptions.None
        };

        BuildReport report = BuildPipeline.BuildPlayer(options);
        BuildSummary summary = report.summary;

        if (summary.result == BuildResult.Succeeded)
        {
            WriteApacheHtaccess(OutputPath);
            Debug.Log($"[WebGLBuilder] Build succeeded: {summary.totalSize} bytes -> {OutputPath}");
            if (Application.isBatchMode)
            {
                EditorApplication.Exit(0);
            }
        }
        else
        {
            Debug.LogError($"[WebGLBuilder] Build failed: {summary.result} (errors: {summary.totalErrors})");
            if (Application.isBatchMode)
            {
                EditorApplication.Exit(1);
            }
        }
    }

    private static void WriteApacheHtaccess(string outputDir)
    {
        string htaccessPath = Path.Combine(outputDir, ".htaccess");
        const string content =
@"<IfModule mod_mime.c>
  AddType application/wasm .wasm
  AddType application/javascript .js
  AddType application/octet-stream .unityweb
</IfModule>

<IfModule mod_dir.c>
  DirectorySlash On
</IfModule>
";
        File.WriteAllText(htaccessPath, content);
        Debug.Log($"[WebGLBuilder] Wrote Apache .htaccess -> {htaccessPath}");
    }
}
