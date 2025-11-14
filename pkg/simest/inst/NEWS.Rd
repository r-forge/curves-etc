% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/simest/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{simest}}% MM: look into ../svn-log-from.all

\section{Changes in version 0.4-1 (2025-11-08)}{
  \subsection{New Features}{
    \itemize{
      \item Improved flexibility \emph{and} informativeness for several
      \code{plot()} methods.

      \item The model fitting functions are no longer S3 generic with a
      \code{default} method, but regular functions (more efficiently).
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item in \code{sim.est()} and \code{simestgcv()}, \code{A[,-c(1,2)]}
      needed to get a \code{ , drop=FALSE } in order to work correctly for
      e.g., \code{n = 6}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Added \pkg{simest} to R-forge, notably \I{Martin Maechler}'s
      \file{curves-etc} project; also adding a \code{BugReports: } entry to
      \file{DESCRIPTION}.

      \item Using \code{match.arg(methods)} and \code{switch(method, ..)} is
	more efficient and more readable code.

      \item added relevant \code{example(fastmerge)}, showing merging with
      default and more merging with larger \code{tol}.

      \item fixed several \file{man/*.Rd} \samp{\\link[pkg]{..}} links,
      also other small \file{*.Rd} issues.

      \item added \file{./tests/} directory and \dQuote{regression} tests,
      which found small inconsistencies, e.g., some methods needing more
      than \eqn{n = 2} observations.
    }
  }
}

%% ---------------------- \pkg{simest} as was archived on CRAN --------
\section{Earlier History and CRAN Archiving of version 0.4}{
  \subsection{CRAN \pkg{simest} Archive}{
    \itemize{
      \item CRAN's archive for \CRANpkg{simest} at
      \url{https://cran.r-project.org/src/contrib/Archive/simest/} shows \preformatted{
	Name        Last modified       Size
	simest_0.1.tar.gz     2015-09-20 22:35        15K
	simest_0.2.tar.gz     2015-10-19 18:17        15K
	simest_0.3.tar.gz     2016-09-05 16:15        22K
	simest_0.4.tar.gz     2017-04-25 16:35        23K}
    }
  }
  \subsection{Archiving Text}{
    \itemize{
      \item Package \pkg{simest} was removed from the CRAN repository.

      \item Formerly available versions can be obtained from the
      \href{https://cran.r-project.org/src/contrib/Archive/simest/}{archive}.

      \item Archived on 2025-03-24 as email to the maintainer is undeliverable.

      \item A summary of the most recent check results can be obtained from the
      \href{https://cran-archive.r-project.org/web/checks/2025/2025-03-24_check_results_simest.html}{check results archive}.

      \item Please use the canonical form
      \url{https://CRAN.R-project.org/package=simest} to link to this page.
    }
  }
}


\section{Changes in version 0.4-0 (2017-04-25)}{
  \subsection{New Features}{
    \itemize{
      \item \code{sim.est()} gets \code{force = TRUE} option and then uses
      the new \code{cvx.lse.con.reg()} which calls \code{conreg()} from package
      \CRANpkg{cobs}.  Hence \pkg{cobs} is now in \pkg{simest}'s reverse
      dependencies.

      \item Now using faster \code{.C(<name>, ....)} instead of
      \code{.C("<name>", ...., PACKAGE = "simest")}
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Documentation: changing \eqn{(x, y)} to \eqn{(t, z)} where appropriate.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item using new \file{src/init.c} and \file{src/simest.h} files
      (necessary for the \code{.C(<name>, ...)} interface.

      \item The C routines called by \code{.C(<name>, *)} documented now
      (even though that is a bit unusual).
    }
  }
}

\section{Changes in version 0.3 (2016-09-05)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
  \subsection{Misc}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.2 (2015-10-19)}{
  \subsection{New Features}{
    \itemize{
      \item .
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item .
    }
  }
  \subsection{Misc}{
    \itemize{
      \item .
    }
  }
}

\section{First version 0.1 (2015-09-20)}{
  \subsection{New Features}{
    \itemize{
      \item .
      \item .
      \item .
      \item .
      \item .
    }
  }
}


