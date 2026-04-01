import codecs
import re

def markdown_to_latex():
    with codecs.open("閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）.md", "r", "utf-8") as f:
        md_text = f.read()

    tex_header = r"""\documentclass[11pt]{ltjsarticle}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{geometry}
\geometry{a4paper, margin=2.5cm}
\usepackage{listings}
\usepackage{xcolor}

\begin{document}

\title{閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）}
\author{木原範昭}
\date{2026年4月}

\maketitle

"""
    tex_footer = "\n\\end{document}\n"

    lines = md_text.split('\n')
    out_lines = []
    
    in_code = False
    in_list = False
    in_math_block = False

    for i in range(len(lines)):
        line = lines[i].rstrip()
        
        if line.startswith('```python'):
            in_code = True
            out_lines.append(r'\begin{lstlisting}[language=Python, basicstyle=\ttfamily\small, breaklines=true]')
            continue
        if line.startswith('```') and in_code:
            in_code = False
            out_lines.append(r'\end{lstlisting}')
            continue
            
        if in_code:
            out_lines.append(line)
            continue
            
        # Headers
        if line.startswith('# '):
            continue # Skip main title, already in \title
        elif line.startswith('## '):
            title = line[3:].strip()
            # remove numbers like "1. " from title
            title = re.sub(r'^\d+\.\s*', '', title)
            out_lines.append(r'\section{' + title + '}')
            continue
        elif line.startswith('### '):
            title = line[4:].strip()
            title = re.sub(r'^\d+\.\d+\s*', '', title)
            out_lines.append(r'\subsection{' + title + '}')
            continue

        # Lists (simple dash)
        is_list_item = line.startswith('- ')
        if is_list_item and not in_list:
            in_list = True
            out_lines.append(r'\begin{itemize}')
        elif not is_list_item and in_list and line.strip() == '':
            # check if next line is list or if list ends
            pass # delay closing
        elif not is_list_item and in_list and line.strip() != '':
            in_list = False
            out_lines.append(r'\end{itemize}')

        if is_list_item:
            content = line[2:]
            out_lines.append(r'\item ' + content)
            continue
            
        # Math blocks
        if line.strip() == '$$':
            if in_math_block:
                out_lines.append(r'\]')
                in_math_block = False
            else:
                out_lines.append(r'\[')
                in_math_block = True
            continue
            
        if line.startswith('$$') and line.endswith('$$') and len(line) > 4:
            out_lines.append(r'\[' + line[2:-2] + r'\]')
            continue

        # Images
        img_match = re.match(r'!\[.*?\]\((.*?)\)', line)
        if img_match:
            img_file = img_match.group(1)
            # Peek next line for caption
            caption = ""
            if i+1 < len(lines) and lines[i+1].startswith('*図'):
                caption_line = lines[i+1]
                caption_match = re.search(r'\*図(\d+):\s*(.*?)\*$', caption_line)
                if caption_match:
                    caption = caption_match.group(2)
            
            out_lines.append(r'\begin{figure}[htbp]')
            out_lines.append(r'  \centering')
            out_lines.append(r'  \includegraphics[width=0.8\textwidth]{' + img_file + '}')
            if caption:
                out_lines.append(r'  \caption{' + caption.replace('_', r'\_') + '}')
            out_lines.append(r'\end{figure}')
            continue
            
        # Skip caption line since handled above
        if line.startswith('*図'):
            continue
            
        # Bold
        line = re.sub(r'\*\*(.*?)\*\*', r'\\textbf{\1}', line)

        out_lines.append(line)

    if in_list:
        out_lines.append(r'\end{itemize}')

    tex_content = tex_header + '\n'.join(out_lines) + tex_footer
    
    with codecs.open("閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）.tex", "w", "utf-8") as f:
        f.write(tex_content)
    print("TeX file created successfully.")

if __name__ == '__main__':
    markdown_to_latex()
