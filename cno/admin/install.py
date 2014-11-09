from biokit.rtools import RPackageManager


__all__ = ["install_all_cellnopt_dependencies"]


def install_all_cellnopt_dependencies(verbose=True):
    """script to install all dependencies for CellNOptR packages

        >>> from cno.admin import install_all_cellnopt_dependencies
        >>> install_all_cellnopt_dependencies()

    """
    pm = RPackageManager()
    installed_packages = pm.installed.index

    # requires those packages to be installed !!
    #packages = ['CellNOptR', 'CNORdt', 'CNORode', 'CNORfeeder', 'CNORfuzz']
    #to_install = []
    #for package in packages:
    #    depends = pm.packages['Depends'].ix[package]
    #    suggests = pm.packages['Suggest'].ix[package]


    cellnopt = ["hash", "Rgraphviz", "RBGL", "graph", "RUnit", "igraph", "XML", "ggplot2", "RCurl"]
    ode= ["Rsolnp", "snowfall", "genalg"]
    feeder = ["catnet","minet"]
    meigor = ["Rsge"]
    fuzzy = ["xtable", "nloptr"]
    dt = []

    packages = sorted(list(set(cellnopt + meigor + fuzzy + feeder + ode + dt)))

    for package in packages:
        if package not in installed_packages:
            if self.verbose:
                print("Installing %s " % package)
            pm.install(package)
        else:
            if self.verbose:
                print("%s already installed. skipped" % package)

    # Rsge not maintained anymore so need to get it from arhive
    if "Rsge" not in installed_packages:
        pm.install_packages("http://cran.r-project.org/src/contrib/Archive/Rsge/Rsge_0.6.3.tar.gz")
    else:
        if self.verbose:
            print("%s already installed. skipped" % "Rsge")

if __name__ == "__main__":
    install_all_cellnopt_dependencies()
