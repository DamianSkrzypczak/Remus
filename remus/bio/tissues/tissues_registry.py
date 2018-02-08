import os


class TissuesFilesRegistry:
    def __init__(self, tissues_dir=os.path.join("db", "tissues"), extension=".bed"):
        self._tissues_dir = tissues_dir
        self._available_tissues = [t for t in os.listdir(tissues_dir) if t.endswith(extension)]
        self._extension = extension

    @property
    def available_tissues(self):
        # TODO: docstring
        return [filename.split(self._extension)[0] for filename in self._available_tissues if
                filename.endswith(self._extension)]

    def get_tissue_path(self, name):
        # TODO: docstring
        if name in self.available_tissues:
            return os.path.join(self._tissues_dir, "{}{}".format(name, self._extension))
