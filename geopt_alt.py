def typename(obj):
    """Returns the type of obj as a string. More descriptive and specific than
type(obj), and safe for any object, unlike __class__."""
    if hasattr(obj, '__class__'):
        return getattr(obj, '__class__').__name__
    else:
        return type(obj).__name__


class GeoPt(object):
  """A geographical point, specified by floating-point latitude and longitude
  coordinates. Often used to integrate with mapping sites like Google Maps.
  May also be used as ICBM coordinates.

  This is the georss:point element. In XML output, the coordinates are
  provided as the lat and lon attributes. See: http://georss.org/

  Serializes to '<lat>,<lon>'. Raises BadValueError if it's passed an invalid
  serialized string, or if lat and lon are not valid floating points in the
  ranges [-90, 90] and [-180, 180], respectively.
  """
  lat = None
  lon = None

  def __init__(self, lat, lon=None):
    if lon is None:

      try:
        split = lat.split(',')
        lat, lon = split
      except (AttributeError, ValueError):
        raise ValueError(
          'Expected a "lat,long" formatted string; received %s (a %s).' %
          (lat, typename(lat)))
    #import pdb;pdb.set_trace()
    try:
      lat = float(lat)
      lon = float(lon)
      if abs(lat) > 90:
        raise ValueError(
          'Latitude must be between -90 and 90; received %f' % lat)
      if abs(lon) > 180:
        raise ValueError(
          'Longitude must be between -180 and 180; received %f' % lon)
    except (TypeError, ValueError):

      raise ValueError(
        'Expected floats for lat and long; received %s (a %s) and %s (a %s).' %
        (lat, typename(lat), lon, typename(lon)))

    self.lat = lat
    self.lon = lon

  def _convert_other(self, other):
    if not isinstance(other, GeoPt):
      try:
        return GeoPt(other)
      except ValueError:
        return NotImplemented
    return other

  def __eq__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    return (self.lat, self.lon) == (other.lat, other.lon)

  def __ne__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    return (self.lat, self.lon) != (other.lat, other.lon)

  def __lt__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    if self.lat < other.lat:
      return True
    if self.lat == other.lat and self.lon < other.lon:
      return True
    return False

  def __le__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    return self < other or self == other

  def __gt__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    if self.lat > other.lat:
      return True
    if self.lat == other.lat and self.lon > other.lon:
      return True
    return False

  def __ge__(self, other):
    other = self._convert_other(other)
    if other is NotImplemented:
      return NotImplemented
    return self > other or self == other

  def __hash__(self):
    """Returns an integer hash of this point.

    Implements Python's hash protocol so that GeoPts may be used in sets and
    as dictionary keys.

    Returns:
      int
    """
    return hash((self.lat, self.lon))

  def __repr__(self):
    """Returns an eval()able string representation of this GeoPt.

    The returned string is of the form 'datastore_types.GeoPt([lat], [lon])'.

    Returns:
      string
    """
    return 'datastore_types.GeoPt(%r, %r)' % (self.lat, self.lon)

  def __unicode__(self):
    return u'%s,%s' % (unicode(self.lat), unicode(self.lon))

  __str__ = __unicode__

  def ToXml(self):
    return u'<georss:point>%s %s</georss:point>' % (unicode(self.lat),
                                                    unicode(self.lon))
