
from lipyd import settings
from lipyd import sample
from lipyd import sampleattrs

peaksfile = settings.get('peaks_example')
peaksfile = 'stard10_pos.csv'

reader = sample.SampleReader(
    input_type = 'peaks',
    fname = peaksfile
)

samples = reader.get_sampleset(
    sampleset_args = {
        'sample_id_proc': sampleattrs.plate_sample_id_processor(),
    }
)

# m/z's for one sample:
samples.mzs_by_sample[:,samples.attrs.sample_id_to_index[('A', 10)]]
# intensities for one sample:
samples.intensities[:,samples.attrs.sample_id_to_index[('A', 10)]]
# RTs for one sample:
samples.rts[:,samples.attrs.sample_id_to_index[('A', 10)]]
