import pylab


__all__ = ["Report", "ReportBool", "ReportODE", "ReportFuzzy" ]




class Report(object):

    def __init__(self):
        pass
    def _make_filename(self, filename):
        return self.report_directory + os.sep + filename
    def savefig(self, filename):
        pylab.savefig(self._make_filename(filename))
    def _get_report_directory(self, directory=None):
        if directory==None:
            directory = self.report_directory
        return directory

    def _init_report(self, directory=None):
        """create the report directroy and return the directory name"""
        directory = self._get_report_directory(directory)
        try:
            os.mkdir(directory)
        except Exception:
            txt = "Existing directory {}. Files may be overwritten".format(self.report_directory)

        if self._overwrite_report == True:
            self.warning(txt)
        else:
            raise IOError(txt)

        for filename in ["dana.css", "reset.css", "tools.js"]:
            filename = gsf("cellnopt.pipeline", "data", filename)
            shutil.copy(filename, directory)
        return directory

    def get_html_reproduce(self):
        text = """
        <p>In a shell, go to the report directory (where is contained this report)
        and either execute the file called rerun.py or typoe these commands in 
        a python shell 
        </p>

        <pre>
        from cellnopt.pipeline import *
        c = CNObool(config=config.ini)
        c.gaBinaryT1()
        c.report()
        </pre>

        <p>You will need the configuration file that was used in the report 
        (<a href="config.ini">config.ini</a>) </p>

        <p>New results will be put in a sub directory so that your current 
        report is not overwritten</p>
        """
        return text

    def plot_fitness(self, show=True, save=False, loglog=True):
        # Requires data Best_score to be found in allresults['fitness']
        pylab.clf()
        self._plot_fitness()
        pylab.ylabel("Fitness")
        if save:
            self.savefig("fit_over_time.png")
            if loglog:
                pylab.semilogy()
                if save:
                    self.savefig("fit_over_time_logy.png")
                pylab.semilogx()
                if save:
                    self.savefig("fit_over_time_logx.png")
        if show == False:
            pylab.close()


class ReportBool(Report):
    def __init__(self):
        super(ReportBool, self).__init__()


class ReportODE(Report):
    def __init__(self):
        super(ReportODe, self).__init__()


class ReportFuzzy(Report):
    def __init__(self):
        super(ReportFuzzy, self).__init__()






