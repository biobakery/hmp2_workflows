def fname_search(collection, key):
    def _iter():
        for item in collection:
            if key in item:
                yield item
        yield ""
    return next(_iter())

