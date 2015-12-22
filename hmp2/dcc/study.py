import cutlass

from .project import default_mixs_dict

class settings:
    name = "ibdmdb" # TODO: remove test when ready to go live
    description = """A new three-year "multi-omic" study investigating the roles played by microbes and their interactions with the human body."""
    contact = "Randall Schwager <schwager@hsph.harvard.edu>"
    center = "Broad Institute"
    subtype = "ibd"
    mixs = dict( biome = "Human, gut",
                 project_name = ("Inflammatory Bowel Disease Multi-omics"
                                 " Database (IBDMDB)") )
    

def load(session):
    res = session.get_osdf().oql_query(
        "ihmp", '"study"[node_type]'
    )['results']
    for data in res:
        if data['meta']['name'] == settings.name:
            return cutlass.Study.load_study(data)
    raise Exception("Unable to find project "+settings.name)


def sync(study):
    study.name = settings.name
    study.description = settings.description
    study.contact = settings.contact
    study.center = settings.center
    study.subtype = settings.subtype
    study.mixs = d = default_mixs_dict()
    d.update(settings.mixs)
    return study

def default(session, project):
    try:
        s = load(session)
    except Exception as e:
        print e
        s = cutlass.Study()
        s.links['part_of'] = [project.id]
    return sync(s)


