import pytest
import os.path as osp
import matplotlib.pyplot as plt
import mpet.plot.plot_data as plot_data

plots = list(plot_data.plotTypes.keys())
plots.remove("text")


@pytest.mark.parametrize("plot", plots)
def test_plot(Dirs, plot):
    refDir, testDir = Dirs
    try:
        plot_data.show_data(osp.join(testDir, "sim_output"), plot, False, False, False)
    except NotImplementedError:
        pass
    except Exception:
        assert False
    plt.close()
    assert True
