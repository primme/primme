set(HEADERS
  primme_eigs_f77.h
  primme_eigs_f90.inc
  primme_eigs.h
  primme_f77.h
  primme_f90.inc
  primme.h
  primme_svds_f77.h
  primme_svds_f90.inc
  primme_svds.h
)

install(FILES ${HEADERS} DESTINATION ${PRIMME_INSTALL_INCLUDEDIR})
