import sys
import os

filename = '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）.md'

with open(filename, 'r', encoding='utf-8') as f:
    content = f.read()

# I corrupted the file with 61041, let's reverse it to $$
resolved_content = content.replace("61041", "$$")

with open(filename, 'w', encoding='utf-8') as f:
    f.write(resolved_content)

print("Restored brackets correctly to $$")
