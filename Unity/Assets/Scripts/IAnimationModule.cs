using TMPro;
using UnityEngine;

public interface IAnimationModule
{
    string Id { get; }
    string Title { get; }
    void Activate(RectTransform panelRoot, Camera mainCamera, TMP_FontAsset font);
    void Deactivate();
}
