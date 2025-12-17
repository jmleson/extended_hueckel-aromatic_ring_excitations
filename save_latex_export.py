
import sympy as sp

from info_document.calculate_with_timeout import calculate_with_timeout
from wrappers import safe_simplify


def save_latex_export(expr: sp.Expr, replacement_text:str="too long to print"):
    expr = try_simplifying(expr=expr)
    text = str(sp.latex(expr, fold_frac_powers=True))
    if len(text) >= 1000:
            text = "\n%" + break_string(text.replace("\n", "\n%")) + "\n"
            text += replacement_text + "\n"
    return text

def break_string(s: str) -> str:
    n= 500
    w = 10
    out, i = [], 0
    while i < len(s):
        t = min(i + n, len(s))
        win = s[max(i, t-w):min(len(s), t+w)]
        p = win.find(" ")
        j = (t-w+p) if p != -1 else t
        out.append(s[i:j])
        i = j + 1
    return "\n%".join(out)


def try_simplifying(expr: sp.Expr):
    simplified_expression = calculate_with_timeout(safe_simplify, (expr,),
                                                   timeout_in_s=120, msg="simplifying")
    if simplified_expression is None:
        return expr
    return simplified_expression