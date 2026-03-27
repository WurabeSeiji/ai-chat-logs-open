import re

def fix_math(text):
    # e^{i\pi S}
    text = re.sub(r'e\s*\n\s*iπS\s*\n', r'e^{i\\pi S}', text)
    text = re.sub(r'e\s*\n\s*iπ/2\s*\n', r'e^{i\\pi/2}', text)
    
    # \neq
    text = re.sub(r'\s*\n\s*=\s*', r' \\neq ', text)
    
    # \sqrt{C}
    text = re.sub(r'平方根（\s*\n\s*C\s*\n\s*\n\s*​\s*\n\s*）', r'$\\sqrt{C}$', text)
    text = re.sub(r'平方根（\s*\n\s*C\s*\n\s*​\s*\n\s*）', r'$\\sqrt{C}$', text)
    
    # matrix P
    text = re.sub(r'P=\(\s*\n\s*0\s*\n\s*E\s*\n\s*E\s*\n\s*0\s*\n\s*​?\s*\n\s*\)', r'P=\\begin{pmatrix} 0 & E \\\\ E & 0 \\end{pmatrix}', text)
    
    # generic symbol_{subscript}^2
    text = re.sub(r'([a-zA-Z+])\s*\n\s*([a-zA-Z]{1,2})\s*\n\s*2\s*\n\s*​?\s*\n?', r'\1_{\2}^2 ', text)
    
    # generic symbol^{2}
    text = re.sub(r'([a-zA-Z+])\s*\n\s*2\s*\n', r'\1^2', text)
    
    # generic symbol_{subscript}
    text = re.sub(r'([a-zA-Z])\s*\n\s*([0-9nijk])\s*\n\s*​?\s*\n?', r'\1_\2', text)
    
    # 180^\circ
    text = re.sub(r'180\s*\n\s*∘\s*\n?', r'180^\\circ ', text)

    # Remove lingering zero width space
    text = text.replace('\u200b', '')
    
    return text

if __name__ == '__main__':
    with open("スピンの定式化方法.md", "r", encoding="utf-8") as f:
        original = f.read()

    fixed = fix_math(original)

    with open("スピンの定式化方法.md", "w", encoding="utf-8") as f:
        f.write(fixed)
    print("Fixed math formulas.")
