SDIR = src
IDIR = include
CC = gcc

PROJECT_BINARY_NAME = Wind_songweiwei

PROJECT_SOURCES = $(SDIR)/ert_main.c
PROJECT_SOURCES += $(SDIR)/Wind_songweiwei_data.c
PROJECT_SOURCES += $(SDIR)/Wind_songweiwei.c
PROJECT_SOURCES += $(SDIR)/rt_nonfinite.c
PROJECT_SOURCES += $(SDIR)/rtGetInf.c
PROJECT_SOURCES += $(SDIR)/rtGetNaN.c
PROJECT_SOURCES += $(SDIR)/rt_matrx.c
PROJECT_SOURCES += $(SDIR)/rt_printf.c
PROJECT_SOURCES += $(IDIR)/sfun_discreteVariableDelay.c
PROJECT_SOURCES += $(SDIR)/rt_look.c
PROJECT_SOURCES += $(SDIR)/rt_look1d.c

IDIR = include
CFLAGS = -I.
CFLAGS += -I$(IDIR)

all:	$(PROJECT_BINARY_NAME)

$(PROJECT_BINARY_NAME):	$(PROJECT_SOURCES)
	$(CC) $(CFLAGS) -o $(PROJECT_BINARY_NAME) $(PROJECT_SOURCES)

clean:
	rm -f $(PROJECT_BINARY_NAME)
