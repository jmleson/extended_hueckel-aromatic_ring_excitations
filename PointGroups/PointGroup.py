

class PointGroup:

    def __init__(self, order:int, n:int):
        self.order = order
        self.n = n
        self.operations = [] # symmetry operations
        self.irreducible_representations = []

        self.total_symmetric_representation = "A1"
        self.dipole_operator_symmetry = None


    def group_order(self):
        return sum([op.amount for op in self.operations])

    def print_multiplication_table(self):
        possible_irreducible_representations = [i.name for i in self.irreducible_representations]
        n = len(possible_irreducible_representations)
        table = [[""] + possible_irreducible_representations]
        for irred_1 in possible_irreducible_representations:
            row = [irred_1]
            for irred_2 in possible_irreducible_representations:
                try:
                    result = self.multiply(irred_1, irred_2)
                    result_str = " + ".join(result)
                except Exception:
                    result_str = "–"
                row.append(result_str)
            table.append(row)
        self.print_table(table)

    def set_up_irreducible_representations(self):
        raise Exception("to-be-implemented in subclass")

    def set_up_symmetry_operations(self):
        raise Exception("to-be-implemented in subclass")

    def multiply(self, irred_1: str, irred_2: str) -> list[str]:
        raise Exception("to-be-implemented in subclass")

    def print_table(self, table):
        BOLD = "\033[1m"
        END = "\033[0m"
        n_cols = len(table[0])
        col_widths = [0] * n_cols
        for row in table:
            for c, cell in enumerate(row):
                col_widths[c] = max(col_widths[c], len(str(cell)))  # nur sichtbare Länge nötig

        for r, row in enumerate(table):
            formatted_cells = []
            for c, cell in enumerate(row):
                text = str(cell)
                text_filled = f"{text:<{col_widths[c]}}"
                if r == 0 or c == 0:
                    text_filled = f"{BOLD}{text_filled}{END}"
                formatted_cells.append(text_filled)

            formatted_row = " │ ".join(formatted_cells)
            print(formatted_row)

            if r == 0:
                total_len = sum(col_widths) + 3 * (n_cols - 1)  # 3 = Länge " │ "
                print("─" * total_len)





