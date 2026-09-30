# EM 自己無撞着有限履歴モデル v1

目的: dipole 放射近似を使わず、二点電荷の retarded Lienard-Wiechert 相互場と Lorentz-Dirac (LAD) 自己力を有限履歴窓から同時に解くための参照実装。

主方程式:
- Maxwell の retarded point-charge solution (Lienard-Wiechert fields)
- Lorentz force
- Lorentz-Dirac radiation reaction (covariant point-charge self-force; 数値実装は3-vector jerk form)

注意:
- LAD は点電荷極限の標準的な自己力方程式だが runaway / preacceleration を持ち得る。
- 本コードは『閉形式解析解』ではなく、厳密方程式を有限履歴上で自己無撞着に反復する数値基盤。
- 重力は含めない。
- harmonic history layer は、履歴窓を波内部位相セルとして扱う将来拡張の接口だけを用意する。
