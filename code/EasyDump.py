
class EasyDump:
    def __init__(self, fn, dec_places):
        """
        Dump variables to a file
        :param fn:
        :param dec_places:
        """
        self.fn = fn
        self.dec_places = dec_places
        self.ff = '{:.' + str(dec_places) + '}'
        self.f = None


    def open(self):
        """
        Open the file
        :return:
        """
        self.f = open(self.fn, "wt")


    def close(self):
        """
        Close the file
        :return:
        """
        self.f.close()


    def write(self, num, name):
        """
        Write to the file
        :param num:         number to write to the file
        :param name:        name of the variable to write
        :return:
        """
        out = name + ':' + self.ff.format(num) + '\n'
        self.f.write(out)


