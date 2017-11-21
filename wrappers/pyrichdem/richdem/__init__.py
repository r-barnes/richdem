import _richdem

def fillDepressions(
  dem,
):
  """Fills all depressions in an elevation model.

     Parameters:
     dem -- A NumPy elevation model

     Returns:
     Modified dem in-place to remove depressions
  """
  return _richdem.rdFillDepressions(dem)

