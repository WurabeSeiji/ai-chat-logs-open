import re
import sys

def format_math(text):
    # 1. Fix block fractions and complex multiline MathJax artifacts
    text = re.sub(r'N=\s*\n\s*4l\s*\n\s*p\s*\n\s*2\s*\n*​\s*\n*\s*A\s*\n*​', r'$$ N = \\frac{A}{4l_p^2} $$', text)
    
    # 2. Fix fragmented superscripts (e.g., L \n 2 \n -> $L^2$)
    # We carefully target typical variables and numbers.
    text = re.sub(r'\b([A-Za-z0-9]+)\s*\n+\s*([0-9]+)\s*\n+(?=[^\s])', r'$\1^{\2}$ ', text)
    text = re.sub(r'\b([A-Za-z0-9]+)\s*\n+\s*([0-9]+)\s*\n+', r'$\1^{\2}$ ', text)
    
    # Also fix sub/superscripts that might have left hanging equals, e.g., $L^2$ = 0 -> $L^2 = 0$
    text = re.sub(r'\$([A-Za-z0-9]+)\^\{([0-9]+)\}\$\s*=\s*(.*?)(?=\s|。|、)', r'$\1^{\2} = \3$', text)
    
    # Fix specific stray pieces where newline broke the equation
    text = re.sub(r'([A-Za-z0-9]+)\s*\n+\s*=\s*([0-9]+)', r'$\1 = \2$', text)
    text = re.sub(r'=\s*\n+\s*([gxy0-9])', r'= \1', text)

    # 3. Format key repeating inline equations
    eqs_to_format = [
        (r'L\s*=\s*X\s*\+\s*Y\s*\+\s*Z\s*-\s*T\s*\+\s*W', r'$L = X + Y + Z - T + W$'),
        (r'L\s*=\s*X\s*\+\s*Y\s*\+\s*Z\s*−\s*T\s*\+\s*W', r'$L = X + Y + Z - T + W$'),
        (r'L\s*=\s*X\s*\+\s*Y\s*\+\s*Z\s*-\s*T\s*\+\s*R', r'$L = X + Y + Z - T + R$'),
        (r'n\s*=\s*T\s*-\s*W', r'$n = T - W$'),
        (r'n\s*=\s*T\s*−\s*W', r'$n = T - W$'),
        (r'0\s*=\s*n\s*-\s*T\s*\+\s*W', r'$0 = n - T + W$'),
        (r'0\s*=\s*X\^2\s*\+\s*Y\^2\s*\+\s*Z\^2\s*-\s*T\^2', r'$0 = X^2 + Y^2 + Z^2 - T^2$'),
        (r'L\^2\s*=\s*X\^2\s*\+\s*Y\^2\s*\+\s*Z\^2\s*-\s*T\^2', r'$L^2 = X^2 + Y^2 + Z^2 - T^2$'),
        (r'-\(電子の質量\)\^2\s*=\s*X\^2\s*\+\s*Y\^2\s*\+\s*Z\^2\s*-\s*T\^2', r'$-(電子の質量)^2 = X^2 + Y^2 + Z^2 - T^2$')
    ]
    for pattern, replacement in eqs_to_format:
        # Avoid replacing already formatted ones by checking if they don't have $
        text = re.sub(r'(?<!\$)\b' + pattern + r'\b(?!\$)', replacement, text)
        text = re.sub(r'(?<!\$)' + pattern + r'(?!\$)', replacement, text) # Without word boundaries for symbols
        
    # Deduplicate $$
    text = text.replace('$$L = X + Y + Z - T + W$', '$L = X + Y + Z - T + W$')
    text = text.replace('$0 = X^2 + Y^2 + Z^2 - T^2$', '$$ 0 = X^2 + Y^2 + Z^2 - T^2 $$')
    text = text.replace('$L^2 = X^2 + Y^2 + Z^2 - T^2$', '$$ L^2 = X^2 + Y^2 + Z^2 - T^2 $$')
    text = text.replace('$-(電子の質量)^2 = X^2 + Y^2 + Z^2 - T^2$', '$$ -(電子の質量)^2 = X^2 + Y^2 + Z^2 - T^2 $$')
    text = text.replace('$$$L', '$$ L')
    text = text.replace('$$$', '$$')

    # Remove instances of spaces between $ and math incorrectly
    text = re.sub(r'\$\s*\$([^$]+)\$\s*\$', r'$$ \1 $$', text)

    return text

if __name__ == "__main__":
    filepath = '波長の次元解析.md'
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    new_content = format_math(content)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(new_content)
        
    print("Formatting complete.")
