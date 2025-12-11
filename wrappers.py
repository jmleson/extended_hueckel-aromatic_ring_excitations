

# Wrapper-Funktionen
def safe_simplify(expr):
    return expr.simplify(doit=False)

def safe_ratsimp(expr):
    return expr.ratsimp()

def safe_expand(expr):
    return expr.expand()