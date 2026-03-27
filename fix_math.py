import re
import sys

def inspect_file(filename):
    with open(filename, "r", encoding="utf-8") as f:
        text = f.read()

    lines = text.split('\n')
    for i in range(1, len(lines)-1):
        prev_line = lines[i-1].rstrip()
        curr = lines[i].strip()
        
        # zero-width spaces
        curr_clean = curr.replace('\u200B', '').strip()
        
        if prev_line and not prev_line.endswith('。') and not prev_line.endswith('、'):
            if curr_clean in ['2', 'n', 'xy', 'yz', 'iπ/2', 'iπS', '∘', '', 'C'] or curr_clean.isdigit():
                print(f"Match around line {i+1}:")
                print(repr(prev_line))
                print(repr(curr))
                print(repr(lines[i+1].strip()))
                print("-" * 20)

if __name__ == "__main__":
    inspect_file("スピンの定式化方法.md")
