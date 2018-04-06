Summary: Program to analyze Whipple Gamma-ray Telescope data files
Name: wuparam
Version: 1.9
Packager: Karl Kosack (kosack@hbar.wustl.edu)
Release: 2
URL: http://jelley.wustl.edu
Source0: %{name}-%{version}.tar.gz
License: GPL
Group: Applications/Science
BuildRoot: %{_tmppath}/%{name}-root
BuildRequires: HDF >= 4.1r2-4
BuildRequires: gsl-devel >= 1.1
BuildRequires: plotutils >= 2.4.1
BuildRequires: XFree86-devel
Requires: plotutils >= 2.4.1-2
Requires: HDF >= 4.1r2-4
Requires: gsl >= 1.1
Requires: XFree86-libs

%description
This is the Washington University Whipple analysis package. Use
it to analyze raw datafiles and generate useful output.
%prep
%setup -q

%build
%configure
make
%install
rm -rf $RPM_BUILD_ROOT
%{makeinstall}

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc AUTHORS NEWS README ChangeLog INSTALL COPYING sample.wuplot wu2txt_data_format.txt
%{_bindir}/wuparam
%{_bindir}/wucut
%{_bindir}/wufit
%{_bindir}/wuopt
%{_bindir}/wu2txt
%{_bindir}/wuplot
%{_bindir}/wuarray
%{_datadir}/wuparam/490.camera
%{_datadir}/wuparam/331.camera
%{_datadir}/wuparam/151.camera
%{_datadir}/wuparam/109.camera
%{_datadir}/wuparam/SuperCuts2001.conf
%{_datadir}/wuparam/SuperCuts2000.conf
%{_datadir}/wuparam/SuperCuts95.conf
%{_datadir}/wuparam/ZCuts2001.conf
%{_datadir}/wuparam/EZCuts2004.conf
%{_datadir}/wuparam/EZCutsForSpectrum.conf
%{_datadir}/wuparam/NoCuts.conf
%{_datadir}/wuparam/SpectralCuts.conf
%{_datadir}/wuparam/alpha-total.gpl
%{_datadir}/wuparam/gamma-rate.gpl
%{_datadir}/wuparam/wuparam-plots.gpl
%{_datadir}/wuparam/image2d.pro
%{_datadir}/wuparam/SAO.edb
%{_datadir}/wuparam/YBS.edb
%{_datadir}/wuparam/animate-by-run.pl
%{_datadir}/wuparam/wuconf.pl
%{_mandir}/man1/wuparam.1.gz
%{_mandir}/man1/wuplot.1.gz
%attr(0755,root,root) %{_datadir}/wuparam/gen-cut-plots.sh
%attr(0755,root,root) %{_datadir}/wuparam/gen-diag-plots.sh


%changelog
* Wed Aug 10 2005 Karl Kosack <kosack@hbar.wustl.edu>
- Updated to include wuarray junk. Also removed wupeds/wugains.

* Sun Jan 31 2005 Karl Kosack <kosack@hbar.wustl.edu>
- Updated spec to work with FC3

* Sun Jan 12 2003 Karl Kosack <kosack@localhost.localdomain>
- See the ChangeLog file for detailed info

* Thu Aug 29 2002 Karl Kosack <kosack@hbar.wustl.edu>
- Version 0.9, now mostly works

* Wed May 15 2002 Karl Kosack <kosack@hbar.wustl.edu>
- Initial build. Very preliminary.






