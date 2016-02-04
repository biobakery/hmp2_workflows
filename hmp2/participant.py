import hashlib
import itertools
from os.path import basename, abspath

from .dcc.subject import fields
from .dcc.sixteen import _getter

def participant_hash(s):
    m = hashlib.md5()
    m.update(s)
    return m.hexdigest()[:8]

def create(output_fname, six_args, wgs_args):
    def _records():
        for packed in (six_args, wgs_args):
            meta_fname, otu_fnames, subj_idx, fname_idx, visit_idx = packed
            subj_idx = _getter(subj_idx)
            fname_idx = _getter(fname_idx)
            visit_idx = _getter(visit_idx)
            fnames = dict([(basename(f), abspath(f)) for f in otu_fnames])
            for record in fields(meta_fname):
                subj_no, visit_no = subj_idx(record), visit_idx(record)
                if not visit_no:
                    continue
                else:
                    visit_no = int(visit_no)
                sample_id = "{}.{}".format(subj_no, visit_no)
                hash_id = "{}.{}".format(participant_hash(subj_no), visit_no)
                fname = fnames.get(fname_idx(record))
                if fname:
                    yield ( sample_id, fname, hash_id )


    def _deduped():
        firstnumeric = lambda l: map(int, l[0].split("."))
        lines = sorted(_records(), key=firstnumeric)
        for _, group in itertools.groupby(lines, key=firstnumeric):
            group = list(group)
            if len(group) > 1:
                for line in group:
                    if "wgs" in line:
                        yield line
                        break
            else:
                yield group[0]
                

    with open(output_fname, 'w') as f:
        print >> f, "\t".join(("#SampleID", "DataFile", "Hash"))
        for line in _deduped():
            print >> f, "\t".join(line)


        
                        
                        
                
            
        
    
    
