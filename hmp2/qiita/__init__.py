from .prep import create_16s_prep
from .sample import create_16s_sample

def fname_search(collection, key):
    def _iter():
        for item in collection:
            if key in item:
                yield item
        yield ""
    return next(_iter())

