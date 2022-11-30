Comparison of different models using Dash
=========================================
You can compare the result of different models using the dashboard build with `Dash <https://dash.plotly.com>`_. To create the dashboard, run the ``mpet_plot_app.py`` script (located in the bin folder) and use the ``-d`` argument to provide a directory with model outputs to include the dashbaord. Each model output should be located in a subfolder of the provided directory. For example, to compare the results of all models saved in subfolders of the folder ``history``, run the command:
``bin/mpet_plot_app.py -d history``. It will try to open the dashbaord in your web browser, where the different models can be identified based on their subfolder name.
Running this script requires the following packages to be installed: `dash <https://pypi.org/project/dash/>`_ and `dash_bootstrap_components <https://pypi.org/project/dash-bootstrap-components/>`_.
