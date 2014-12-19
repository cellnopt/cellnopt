# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import os
import shutil


import pylab
import easydev
import pandas as pd
from biokit.rtools import RPackageManager


__all__ = ["Report", "ReportBool", "ReportODE", "ReportFuzzy", "ReportDT"]


class HTMLTable(object):
    def __init__(self, df, name, **kargs):
        self.df = df
        self.name = name
        self.kargs = kargs.copy()

    def to_html(self, index=False):
        pd.set_option('display.max_colwidth', -1)
        kargs = self.kargs.copy()
        kargs['index'] = index
        kargs['bold_rows'] = True
        table = self.df.to_html(**kargs)
        pd.set_option('display.max_colwidth', 50)
        return table


class Report(easydev.Logging):

    def __init__(self, formalism='base', directory='report', tag=None, filename="index.html",
                 overwrite=True, verbose=True, dependencies=True):
        super(Report, self).__init__(verbose)
        self.formalism = formalism
        from cno import version
        self.version = version
        self.tag = tag

        if self.tag is None:
            self.report_directory = "_".join([directory, self.formalism])
        else:
            self.report_directory = "_".join([directory, self.formalism, self.tag])

        self._overwrite_report = overwrite
        self.sections = []
        self.section_names = []
        self.index = filename
        self.Rdependencies = []

    def show(self):
        from browse import browse as bs
        bs(self.report_directory + os.sep + self.index)

    def close_body(self):
        return "</body>"

    def close_html(self):
        return "</html>"

    def get_footer(self):
        return self.close_body() + "\n" + self.close_html()

    def _make_filename(self, filename):
        return self.report_directory + os.sep + filename

    def savefig(self, filename):
        pylab.savefig(self._make_filename(filename))

    def _init_report(self, directory=None):
        """create the report directroy and return the directory name"""
        self.sections = []
        self.section_names = []

        if directory is None:
            directory = self.report_directory

        # if the directory already exists, print a warning
        try:
            os.mkdir(directory)
        except Exception:
            # already exists
            txt = "Existing directory {}. Files may be overwritten".format(self.report_directory)
            if self._overwrite_report is True:
                self.warning("Directory %s exists already. Files may be overwritten" % directory)
            else:
                raise IOError('Directory %s exists already. Set _overwrite_report to True or delete the directory' % directory)

        for filename in ["dana.css", "reset.css", "tools.js"]:
            filename = easydev.gsf("cno", "", filename)
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
     <title>CNO report</title>
     <link rel="stylesheet" href="dana.css" type="text/css" />
     <script type='text/javascript' src='tools.js'></script>
 </head>

 <body>
  <div class="document" id="unset">

     <h1 class="title">CNO%(formalism)s analysis summary</h1>
     <h2 class="subtitle", id="unset2">Report created with cno.%(formalism)s.CNO%(formalism)s (version %(version)s)</h2>
     <p>See <a href="https://www.cellnopt.org">cellnopt homepage</a> for downloads and documentation.</p>""" % params
        return str_

    def get_html_reproduce(self):
        text = """
        <p>In a shell, go to the report directory (where is contained this report)
        and either execute the file called rerun.py or typoe these commands in
        a python shell
        </p>

        <pre>
        from cno import CNObool
        c = CNObool(config=config.ini)
        c.optimise()
        c.onweb()
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

        dependencies = easydev.get_dependencies('cno')
        names = [x.project_name for x in dependencies]
        versions = [x.version for x in dependencies]
        links = ["""https://pypi.python.org/pypi/%s"""%p for p in names]

        if len(self.Rdependencies):
            self._rpm = RPackageManager()
            for package in self.Rdependencies:
                names.append(package)
                versions.append(self._rpm.get_package_version(package))
                links.append("http://www.cellnopt.org")

        df = pd.DataFrame({
            'package': ["""<a href="%s">%s</a>"""%(links[i],p) for i,p in enumerate(names)],
            'version':versions})

        table = HTMLTable(df, name="dependencies", escape=False)

        from cno.misc.dependencies import plot_dependencies
        c = plot_dependencies(show=False, filename=self._make_filename('dependencies.svg'))

        return table


    def add_section(self, content, title, references=[], position=None):

        reftxt = self._create_references(references)
        section = """<div class="section" id="%(id)s">
        <h2> <a class="toc-backref" href="#id%(index)s">%(title)s</a></h2>

        %(references)s\n
        %(content)s
    </div>
        """ % {'title':title, 'references':reftxt,'content':content,
               'id': title.replace(" ", "_"), 'index':len(self.sections)+1}
        # check that it is correct
        if position is not None:
            self.sections.insert(position, section)
            self.section_names.insert(position, title)
        else:
            self.sections.append(section)
            self.section_names.append(title)

    def get_toc(self):
        """
        """
        toc = """<div class="contents local topic" id="contents">
        <ul class="simple">"""
        for i, name in enumerate(self.section_names):
            toc += """<li>
%(i)s - <a class="reference internal" href="%(href)s" id="%(id)s">  %(name)s</a>
</li>""" % {'i':i+1, 'name':name, 'href':"#"+name.replace(" ", "_"), 'id':'id%s' % str(i+1)}
        toc += """</ul>\n</div>"""
        return toc

    def _create_references(self, references):
        if len(references) == 0:
            return ""

        txt = """
        <div class="admonition-documentation admonition">
            <p class="first admonition-title">Documentation</p>
            <ul class="last simple">
            <li>"""

        for ref in references:
            txt += """      <a class="reference external" href=%(url)s>%(title)s</a>""" %  {'url':ref[0],'title':ref[1]}
        txt += """
            </li>
            </ul>
        </div>"""
        return txt

    def write(self, filename):
        fh =  open(self.report_directory + os.sep + filename, "w")

        contents = self.get_header()
        #contents += self.get_toc()

        # Get toc should be done here and no more sections should be added
        self.add_section(self.get_toc(), "Contents", position=0)
        for i, section in enumerate(self.sections):
            if i==0:
                contents += section
            else:
                contents += section.replace("<h2>", "<h2> %s - " %i, 1)

        contents += self.get_table_dependencies().to_html()
        contents += "<hr>" + self.get_time_now()
        contents += self.get_footer()

        import bs4
        contents = bs4.BeautifulSoup(contents).prettify()
        fh.write(contents)
        fh.close()

    def _make_filename(self, filename):
        return self.report_directory + os.sep + filename

    def savefig(self, filename):
        pylab.savefig(self._make_filename(filename))

    def report(self, browse=True):
        self._create_report()
        try:
            self.create_report_images()
        except:
            pass
        if browse:
            self._report.show()

    def _create_report(self):
        raise NotImplementedError

    def _create_report_images(self):
        raise NotImplementedError


class ReportBool(Report):
    def __init__(self):
        super(ReportBool, self).__init__('bool')
        self.Rdependencies.append('CellNOptR')


class ReportODE(Report):
    def __init__(self):
        super(ReportODE, self).__init__('ode')
        self.Rdependencies.append('CNORode')


class ReportFuzzy(Report):
    def __init__(self):
        super(ReportFuzzy, self).__init__('fuzzy')
        self.Rdependencies.append('CNORfuzzy')


class ReportDT(Report):
    def __init__(self):
        super(ReportDT, self).__init__('dt')
        self.Rdependencies.append('CNORdt')


class ReportASP(Report):
    def __init__(self):
        super(ReportASP, self).__init__('sap')


class ReportILP(Report):
    def __init__(self):
        super(ReportILP, self).__init__('ilp')




