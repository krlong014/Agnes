
class XMLHeader:
    def __init__(self, tag):
        self.tag = tag
        self.attributes = {}

    def addAttribute(self, name, val):
        self.attributes[name]=str(val)

    def header(self):
        h = '<' + self.tag
        for k,v in self.attributes.items():
            h = '%s %s=\"%s\"' % (h, k, v)
        h = h + '>'
        return h

    def footer(self):
        return '</%s>' % self.tag


if __name__=='__main__':

    h = XMLHeader('Test')
    h.addAttribute('n', 32)
    h.addAttribute('color', 'blue')
    h.addAttribute('pi', 3.14)

    print(h.header())
    print('bob')
    print(h.footer())
