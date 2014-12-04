import os
import pylab
import easydev
import shutil
from easydev import gsf

__all__ = ["Report", "ReportBool", "ReportODE", "ReportFuzzy" ]




class HTML(object):
    def close_body(self):
        return "</body>"
    def close_html(self):
        return "</html>"
    def get_footer(self):
        return self.close_body() + "\n" + self.close_html()






class HTMLTable(object):
    def __init__(self, df, name):
        self.df = df
        self.name = name
    def open_xml_table(self):
        s = """<Table Name="%s:table">\n""" % self.name
        return s
    def close_xml_table(self):
        s = """</Table>\n"""
        return s
    def write_xml_data(self, d):
        delimiter = self.delimiter + " "
        s = "        " + delimiter.join(['"%s"' % str(x) for x in d]) + "\n"
        return s
    def to_xml(self):
        s = self.open_xml_table()
        s += self.create_xml_columns()
        s += self.open_xml_stream()
        for d in self.data:
            s += self.write_xml_data(d)
        s += self.close_xml_stream()
        s += self.close_xml_table()
        return s
    def create_xml_columns(self):
        s = ""
        for t, c in zip(self._coltypes, self._columns):
            s += """    <Column Name="%s:%s" Type="%s"/>\n""" % (self.name, c, t)
        return s

    def close_xml_stream(self):
        if len(self.data):
            s = """    </Stream>\n"""
        else:
            s = """\n    </Stream>\n"""
        return s

    def open_xml_stream(self):
        s = """    <Stream Name="%s:table" Delimiter="%s">\n""" % (self.name,
                                                                   self.delimiter)
        return s



    def to_html(self):
        return self.df.to_html()




class Report(object):

    def __init__(self, formalism):
        self.formalism = formalism
        from cno import version
        self.version = version

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
            filename = gsf("cno", "data", filename)
            shutil.copy(filename, directory)
        return directory

    def get_header(self):
        """a possible common header ? """
        params = {"formalism": self.formalism,
                   "version": self.version}
        str_ =  """
 <!DOCTYPE html>
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
     <head>
     <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
     <title>CellNOpt pipeline report</title>
     <link rel="stylesheet" href="dana.css" type="text/css" />
     <script type='text/javascript' src='tools.js'></script>
 </head>

 <body>
  <div class="document" id="unset">

     <h1 class="title">CNO%(formalism)s analysis summary</h1>
     <h2 class="subtitle", id="unset2">Report created with cellnopt.pipeline.cno%(formalism)s.CNO%(formalism)s (version %(version)s)</h2>
     <p>See <a href="https://www.cellnopt.org">cellnopt homepage</a> for downloads and documentation.</p>""" % params
        return str_

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


    def get_time_now(self):
        import datetime
        username = os.environ["USERNAME"]
        msg = '<div class="date">Created on ' + str(datetime.datetime.now())
        msg +=  " by " + username +'</div>'
        return msg

    def get_table_dependencies(self):
        """Returns dependencies of the pipeline into a HTML/XML table

        dependencies are the python dependencies as returned by pkg_resource.
        additionally, r dependencies added in :attr:`dependencies` are also added.

        """
        table = HTMLTable(name="dependencies")
        table.add_columns(["software", "revision"], ['str', 'str'])
        for x in easydev.get_dependencies("cellnopt.pipeline"):
            table.add_data([x.project_name, x.version])
        for dep in self.Rdependencies:
            try:
                import rtools
                version = rtools.RPackage(dep).version
                table.add_data([dep, version])
            except:
                table.add_data([dep, "unknown. could not use Rtools/rpy2"])
                #logging.warning("Could not find version of %s" % dep)
        return table


class ReportBool(Report):
    def __init__(self):
        super(ReportBool, self).__init__('bool')

class ReportODE(Report):
    def __init__(self):
        super(ReportODE, self).__init__('ode')


class ReportFuzzy(Report):
    def __init__(self):
        super(ReportFuzzy, self).__init__('fuzzy')


class ReportDT(Report):
    def __init__(self):
        super(ReportDT, self).__init__('dt')


class ReportASP(Report):
    def __init__(self):
        super(ReportASP, self).__init__('sap')


class ReportILP(Report):
    def __init__(self):
        super(ReportILP, self).__init__('ilp')




