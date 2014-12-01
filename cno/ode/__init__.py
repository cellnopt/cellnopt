

def install_dependencies_ode():

    from biokit.rtools import package
    pm = package.RPackageManager()


    if "Rsge" not in pm.installed.index:
        #rtools.install_packages("http://cran.r-project.org/src/contrib/Archive/Rsge/Rsge_0.6.3.tar.gz")
        pm.install("Rsge")
        #pm.install_packages(["snowfall", "Rsolnp"], repos=None)
    if "MEIGOR" not in pm.installed.index:
        pm.install_packages("http://www.cellnopt.org/downloads/MEIGOR_0.99.6_svn3222.tar.gz",
            type="source")


from .cnorode import CNORode
