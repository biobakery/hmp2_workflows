import logging
from itertools import islice
from itertools import chain

import cutlass.iHMPSession

from . import project
from . import study
from . import subject
from . import visit
from . import sixteen
from . import wgs
from . import wms as wts

class settings:
    username = 'rschwager'
    
def submit(hosp_joined_fname, sixt_joined_fname, sixt_fnames,
           wgs_joined_fname, wgs_fnames, wts_joined_fname, wts_fnames,
           dcc_user, dcc_pass):
    logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
    logging.getLogger().setLevel(logging.DEBUG)
    # 2aeafd32a30284042042e3d58787ba42
    session = cutlass.iHMPSession(dcc_user, dcc_pass)
    pr = project.default(session)
    st = study.default(session, pr)
    st.save()
    subjects = list(subject.from_file(hosp_joined_fname, st))
    validation_errors = [ (s, s.validate()) for s in subjects 
                          if not s.is_valid() ]
    if validation_errors:
        raise visit.SaveError(validation_errors)
    has_id = set()
    for s in subjects:
        if not s.id:
            s.save()
        else:
            has_id.add(s)
    for s in has_id:
        s.save()

    visits = map(
        list, visit.from_file(hosp_joined_fname, subjects, nest=True))
    visit.save_all(visits)

    samples = list()
    for grp in chain.from_iterable(visits):
        samples.extend([s for s in islice(grp, 1, None)])

    ps_groups = map(
       list, sixteen.from_file(sixt_joined_fname, samples, sixt_fnames))
    sixteen.save_all(ps_groups)

    ps_groups = map(
       list, sixteen.from_file(
           wgs_joined_fname, samples, wgs_fnames,
           sample_id_idx=lambda r: "{}-{}".format(r[1][1:3], r[1][3:]),
           fname_idx=1,
           parse_record=wgs.parse_record)
    )
    sixteen.save_all(ps_groups)

    ps_groups = map(
       list, sixteen.from_file(
           wts_joined_fname, samples, wts_fnames, sample_id_idx=0, fname_idx=1,
           parse_record=wts.parse_record)
    )
    sixteen.save_all(ps_groups)

    
