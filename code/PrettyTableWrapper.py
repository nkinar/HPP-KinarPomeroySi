from prettytable import PrettyTable
from check_if_iterable import check_if_iterable
from float_round import float_round

class PrettyTableWrapper:

    def __init__(self, dec_places):
        self.dec_places = dec_places
        self.tabd = PrettyTable()

    def field_names(self, x):
        self.tabd.field_names = x

    def add_row(self, x):
        if not check_if_iterable(x):
            return
        row = []
        for k in x:
            if isinstance(k, str):
                row.append(k)
                continue
            else:
                row.append(float_round(k, self.dec_places))
        self.tabd.add_row(row)
    # DONE

    def get_obj(self):
        return self.tabd

    def get_html_string(self):
        return self.tabd.get_html_string()


