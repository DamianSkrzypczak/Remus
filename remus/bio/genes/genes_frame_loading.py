import pandas as pd


# TODO: docstring

class GenesFrameLoader:
    # TODO: docstring
    def __init__(self, source_path, sep="\t"):
        # TODO: docstring
        # TODO: get path and sep from config
        self._source_path = source_path
        self._raw_genes_data_frame = pd.read_csv(self._source_path, dtype=object, sep=sep)
        self._genes_frame = None

    @property
    def genes_frame(self):
        # TODO: docstring
        if self._genes_frame is None:
            self._prepare_genes_frame()
        return self._genes_frame

    def _prepare_genes_frame(self):
        # TODO: docstring
        genes_names = self._extract_names()
        genes_source = self._extract_sources()
        genes_frame = pd.DataFrame({"gene_name": genes_names, "gene_source": genes_source})
        genes_frame.set_index("gene_name", inplace=True)
        genes_frame["gene_name"] = genes_frame.index
        self._genes_frame = genes_frame

    def _extract_sources(self):
        return self._raw_genes_data_frame.iloc[:, [1, 3, 4]].apply(lambda x: '\t'.join(x), axis=1)

    def _extract_names(self):
        return self._raw_genes_data_frame.iloc[:, -1]
