__author__ = 'srio'

import json


class GFile(object):

    def __init__(self,filename=None):

        if filename is None:
            self.gfile_as_dictionary = None
        else:
            self.load_gfile(filename)


    def load_gfile(self,filename,use_brackets_instead_parenthesis=True):
        fp = open(filename) # Open file on read mode
        lines = fp.read().split("\n") # Create a list containing all lines
        fp.close() # Close file

        out_dictionary = {}
        for line in lines:
            mylist = line.split("=")
            if len(mylist) == 2:
                key = mylist[0].strip()
                if use_brackets_instead_parenthesis:
                    key = key.replace("(","[")
                    key = key.replace(")","]")
                value = mylist[1].strip()
                try:
                    if "." in value:
                        out_dictionary[key] = float(value)
                    else:
                        out_dictionary[key] = int(value)
                except:

                    out_dictionary[key] = value

        self.gfile_as_dictionary = out_dictionary


    @classmethod
    def _dictionary2object(cls,dict1):

        class K(object):
            pass

        k = K()
        for key in dict1.keys():
            setattr(k, key, dict1[key])
        return k

    def get_as_dictionary(self):
        return self.gfile_as_dictionary

    def get_as_object(self):
        return self._dictionary2object(self.get_as_dictionary())

    def get_as_json(self):
        return json.dumps(self.get_as_dictionary(), sort_keys=True, indent=4)

if __name__ == "__main__":

    g = GFile()

    g.load_gfile(filename="start.01")



    print("\nDICT: \n",g.get_as_dictionary())

    print("\nOBJECT: \n",dir(g.get_as_object()))

    print(">>>>>Y_ROT",g.get_as_object().Y_ROT )
    print(">>>>>ISTAR1",g.get_as_object().ISTAR1 )

    tmp = getattr(g.get_as_object(),"THICK[6]")
    print(">>>>>THICK[6]",tmp,g.get_as_object())


    print("\nJSON: \n",g.get_as_json())


