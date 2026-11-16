using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.InputSystem;
using UnityEngine.InputSystem.EnhancedTouch;
using Touch = UnityEngine.InputSystem.EnhancedTouch.Touch;

/// <summary>
/// PC（マウスドラッグ・ホイール）とスマホ（1本指スワイプ・2本指ピンチ）の両方に対応した
/// 3D シミュレーション／ゲーム向けオービットカメラコントローラ。
/// </summary>
public class OrbitCameraController : MonoBehaviour
{
    [Header("Target & Distance")]
    public Transform target;
    public Vector3 targetOffset = Vector3.zero;
    public float distance = 8.0f;
    public float minDistance = 2.0f;
    public float maxDistance = 30.0f;

    [Header("Sensitivity")]
    public float mouseOrbitSpeed = 0.25f;
    public float touchOrbitSpeed = 0.2f;
    public float scrollZoomSpeed = 0.01f;
    public float pinchZoomSpeed = 0.02f;

    private float _yaw = 45f;
    private float _pitch = 25f;
    private bool _mouseLockedToUI;
    private bool _touchLockedToUI;

    public void SetView(float yaw, float pitch, float dist, Vector3 focusOffset)
    {
        _yaw = yaw;
        _pitch = Mathf.Clamp(pitch, -80f, 80f);
        distance = Mathf.Clamp(dist, minDistance, maxDistance);
        targetOffset = focusOffset;
        UpdateCameraTransform();
    }

    private void OnEnable()
    {
        EnhancedTouchSupport.Enable();
    }

    private void OnDisable()
    {
        EnhancedTouchSupport.Disable();
    }

    private void Start()
    {
        Vector3 angles = transform.eulerAngles;
        _yaw = angles.y;
        _pitch = angles.x;
        UpdateCameraTransform();
    }

    private void LateUpdate()
    {
        HandleTouchInput();
        HandleMouseInput();
        UpdateCameraTransform();
    }

    private static bool IsPointerOverUI(int pointerId = -1)
    {
        return EventSystem.current != null && EventSystem.current.IsPointerOverGameObject(pointerId);
    }

    private void HandleMouseInput()
    {
        Mouse mouse = Mouse.current;
        if (mouse == null) return;

        if (mouse.leftButton.wasPressedThisFrame)
        {
            _mouseLockedToUI = IsPointerOverUI();
        }
        if (mouse.leftButton.wasReleasedThisFrame)
        {
            _mouseLockedToUI = false;
        }

        if (mouse.leftButton.isPressed && !_mouseLockedToUI)
        {
            Vector2 delta = mouse.delta.ReadValue();
            _yaw += delta.x * mouseOrbitSpeed;
            _pitch = Mathf.Clamp(_pitch - delta.y * mouseOrbitSpeed, -80f, 80f);
        }

        if (!IsPointerOverUI())
        {
            float scroll = mouse.scroll.ReadValue().y;
            if (Mathf.Abs(scroll) > 0.01f)
            {
                distance = Mathf.Clamp(distance - scroll * scrollZoomSpeed, minDistance, maxDistance);
            }
        }
    }

    private void HandleTouchInput()
    {
        var touches = Touch.activeTouches;
        if (touches.Count == 0)
        {
            _touchLockedToUI = false;
            return;
        }

        if (touches.Count == 1)
        {
            var t = touches[0];
            if (t.phase == UnityEngine.InputSystem.TouchPhase.Began)
            {
                _touchLockedToUI = IsPointerOverUI(t.touchId);
            }
            if (!_touchLockedToUI && t.phase == UnityEngine.InputSystem.TouchPhase.Moved)
            {
                _yaw += t.delta.x * touchOrbitSpeed;
                _pitch = Mathf.Clamp(_pitch - t.delta.y * touchOrbitSpeed, -80f, 80f);
            }
        }
        else if (touches.Count >= 2 && !_touchLockedToUI)
        {
            var t0 = touches[0];
            var t1 = touches[1];
            Vector2 prev0 = t0.screenPosition - t0.delta;
            Vector2 prev1 = t1.screenPosition - t1.delta;

            float prevDist = Vector2.Distance(prev0, prev1);
            float currDist = Vector2.Distance(t0.screenPosition, t1.screenPosition);
            float diff = currDist - prevDist;

            distance = Mathf.Clamp(distance - diff * pinchZoomSpeed, minDistance, maxDistance);
        }
    }

    private void UpdateCameraTransform()
    {
        Vector3 focus = (target != null ? target.position : Vector3.zero) + targetOffset;
        Quaternion rot = Quaternion.Euler(_pitch, _yaw, 0f);
        transform.position = focus + rot * new Vector3(0f, 0f, -distance);
        transform.LookAt(focus);
    }
}
