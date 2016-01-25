import cutlass.Project
from cutlass.mixs import MIXS

def default_mixs_dict():
    return dict([ (k, v()) for k, v in MIXS._fields.iteritems() ])

desc = \
"""A new three-year "multi-omic" study investigating the roles played by microbes and their interactions with the human body."""

class settings:
    name = "iHMP"
    
def default(session):
    res = session.get_osdf().oql_query(
        "ihmp", '"project"[node_type]'
    )['results']
    for data in res:
        if data['meta']['name'] == settings.name:
            return cutlass.Project.load_project(data)
    raise Exception("Unable to find project "+settings.name)


