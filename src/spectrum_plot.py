# -- 2019/01/22 --
import os
import importlib

import lipyd.ms2 as ms2
import lipyd.common as common
import lipyd.plot as plot
import lipyd.settings as settings




if __name__ == "__main__":
    importlib.reload(settings)
    importlib.reload(ms2)
    importlib.reload(plot)

    mgfname = 'pos_examples.mgf'
    ionmode = 'pos'
    scan_id = 2397
    scan_args = {}

    mgfpath = os.path.join(common.ROOT, 'data', 'ms2_examples', mgfname)
    scan = ms2.Scan.from_mgf(mgfpath, scan_id, ionmode, **scan_args)

    p = plot.SpectrumPlot(
        mzs = scan.mzs,
        intensities = scan.intensities,
        annotations = scan.annot,
        scan_id = scan.scan_id,
        ionmode = ionmode,
        sample_name = 'BPI_invivo',
        grid_cols = 1
    )

