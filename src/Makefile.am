bin_PROGRAMS = wuparam wucut wufit wuopt wuplot wu2txt wuarray

EXTRA_DIST = Doxyfile wu2root.cpp Makefile.root

AM_LDFLAGS = -L/usr/X11R6/lib 
LDADD = -lmfhdf -ldf -lXaw -lXmu -lXt -lSM -lICE -lXext -lX11 

wupsources = \
	Camera.cpp Camera.h \
	Config.cpp Config.h \
	DataReader.cpp DataReader.h DataRecords.h\
	DataWriter.cpp DataWriter.h \
	GainFinder.cpp GainFinder.h \
	ImageAnalyzer.cpp ImageAnalyzer.h \
	MuonImageAnalyzer.cpp MuonImageAnalyzer.h \
	ImageCleaner.cpp ImageCleaner.h \
	Parameterizer.cpp Parameterizer.h \
	PedestalFinder.cpp PedestalFinder.h \
	ProgressBar.cpp ProgressBar.h \
	Histogram.cpp Histogram.h \
	GaussianDeviate.cpp GaussianDeviate.h \
	Types.h Exceptions.h \
	ParamDataReader.cpp ParamDataReader.h \
	Cutter.cpp Cutter.h \
	StarCatalog.cpp StarCatalog.h \
	Image2D.cpp Image2D.h \
	PlotMaker.cpp PlotMaker.h \
        SuperCutter.cpp SuperCutter.h \
        ZCutter.cpp ZCutter.h \
        EZCutter.cpp EZCutter.h \
        SpectralCutter.cpp SpectralCutter.h \
        EnergySpectrum.cpp EnergySpectrum.h \
	AngleConverters.h \
	ImageFilter.cpp ImageFilter.h \
	Log.cpp Log.h 

wuparam_SOURCES = $(wupsources) wuparam.cpp wuparam.h 
wuparam_LDADD =  $(LDADD)  -lXaw -lXmu -lXt -lSM -lICE -lXext -lX11 
wucut_LDADD =  $(LDADD)  -lXaw -lXmu -lXt -lSM -lICE -lXext -lX11 

#wuparam_mpi_SOURCES = $(wupsources) wuparam-mpi.cpp wuparam.h 
#wuparam_mpi_LDADD = $(wuparam_LDADD) -lmpi -llam

wucut_SOURCES = \
	wucut.cpp \
	Camera.cpp Camera.h \
	Config.cpp Config.h \
	Cutter.cpp Cutter.h \
	DataReader.cpp DataReader.h DataRecords.h \
	DataWriter.cpp DataWriter.h \
	ImageAnalyzer.cpp ImageAnalyzer.h \
	ParamDataReader.cpp ParamDataReader.h \
	ProgressBar.cpp ProgressBar.h \
	StarCatalog.cpp StarCatalog.h \
	Histogram.cpp Histogram.h \
	Types.h Exceptions.h \
	Image2D.cpp Image2D.h \
	EnergySpectrum.cpp EnergySpectrum.h \
        SuperCutter.cpp SuperCutter.h \
        ZCutter.cpp ZCutter.h \
        EZCutter.cpp EZCutter.h \
	PlotMaker.cpp PlotMaker.h \
        SpectralCutter.cpp SpectralCutter.h \
	AngleConverters.h \
	Log.cpp Log.h 

wuarray_SOURCES = \
	$(wupsources) \
	wuarray.cpp

wufit_SOURCES = \
	$(wupsources) \
        wufit.cpp \
        Fitter.cpp Fitter.h 



wuopt_SOURCES =  wuopt.cpp \
	Config.cpp Config.h \
	Cutter.cpp Cutter.h \
	DataReader.cpp DataReader.h DataRecords.h \
	DataWriter.cpp DataWriter.h \
	ImageAnalyzer.cpp ImageAnalyzer.h \
	ParamDataReader.cpp ParamDataReader.h \
	ProgressBar.cpp ProgressBar.h \
	StarCatalog.cpp StarCatalog.h \
	Histogram.cpp Histogram.h \
	Types.h Exceptions.h \
	Image2D.cpp Image2D.h \
	EnergySpectrum.cpp EnergySpectrum.h \
        SuperCutter.cpp SuperCutter.h \
        ZCutter.cpp ZCutter.h \
        EZCutter.cpp EZCutter.h \
	PlotMaker.cpp PlotMaker.h \
	Camera.cpp Camera.h \
        SpectralCutter.cpp SpectralCutter.h \
	Log.cpp Log.h 


wuplot_SOURCES = wuplot.cpp \
	Config.cpp Config.h \
	Image2D.cpp Image2D.h \
	PlotMaker.cpp PlotMaker.h \
	Types.h Exceptions.h \
	AngleConverters.h \
	StarCatalog.cpp StarCatalog.h 

wu2txt_SOURCES = wu2txt.cpp \
	Config.cpp Config.h \
	Types.h Exceptions.h \
	ProgressBar.cpp ProgressBar.h \
	DataReader.cpp DataReader.h DataRecords.h \
	ParamDataReader.cpp ParamDataReader.h \
	ImageAnalyzer.cpp ImageAnalyzer.h \
	Log.cpp Log.h \
	DataWriter.cpp DataWriter.h





