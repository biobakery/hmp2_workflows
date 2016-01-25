import logging
from itertools import islice
from itertools import chain

import cutlass.iHMPSession

from . import project
from . import study
from . import subject
from . import visit
from . import sixteen


class settings:
    username = 'rschwager'
    
def submit(hosp_joined_fname, sixt_joined_fname, sixt_fnames,
           dcc_user, dcc_pass):
    logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")
    logging.getLogger().setLevel(logging.DEBUG)
    # 2aeafd32a30284042042e3d58787ba42
    session = cutlass.iHMPSession(dcc_user, dcc_pass)
    pr = project.default(session)
    st = study.default(session, pr)
    st.save()
    subjects = list(subject.from_file(hosp_joined_fname, st))
    for s in subjects:
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
