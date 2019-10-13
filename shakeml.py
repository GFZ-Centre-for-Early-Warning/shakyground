#Parser for shakemap xml files
import pandas
import lxml.etree as le
import io
import numpy as np
import quakeml
import datetime

#FIXME: shakemap_id event_id

#shakemlfile = 'grid.xml'
#with open(shakemlfile,'r') as f:
#    shakeml = f.read()

class quakemap():
    '''
    Class that keeps all information related to a shakemap
    '''
    def __init__(self,event,event_specific_uncertainty,shakemap,units,regular_grid):
        self.event = event #pandas series
        self.event_specific_uncertainty = event_specific_uncertainty #pandas dataframe
        self.shakemap = shakemap #pandas dataframe
        self.units = units #pandas series
        self.regular_grid = regular_grid#regular or irregular

class site_params():
    '''
    Class that keeps all information related to sites
    '''
    def __init__(self,sites,units,regular_grid):
        self.sites = sites #pandas dataframe
        self.units = units #pandas series
        self.regular_grid = regular_grid #regular or irregular

def shakemap2quakemap(shakemlfile):
    '''
    Reads a shakemal xml file and returns a event table and a shakmap grid
    If no event is provided it assumes that this is a site location file
    '''
    #shakeml = bytes(bytearray(shakeml, encoding='utf-8'))
    #shakeml = le.XML(shakeml)
    #tree = le.parse(shakemlfile)
    #large text node parser
    #read string or file
    try:
        shakeml = le.parse(shakemlfile)
    except:
        #maybe string
        parser = le.XMLParser(huge_tree=True)
        #shakeml = le.parse(io.StringIO(shakemlfile),parser)
        #shakeml = le.parse(shakemlfile,parser)
        try:
            inp = io.BytesIO(shakemlfile)
        except TypeError:
            inp = io.StringIO(shakemlfile)
        shakeml = le.parse(inp,parser)
    nsmap = shakeml.getroot().nsmap
    shakeml = shakeml.getroot()
    #find event
    smevent=shakeml.find('event',namespaces = nsmap)

    #It can be a sites only definition
    siteml=False
    if smevent==None:
        siteml=True
    else:
        #event attributes
        index = [i for i in range(max(1,len(smevent)))]
        columns=['eventID', 'Agency', 'Identifier', 'year', 'month', 'day', 'hour', 'minute', 'second', 'timeError', 'longitude', 'latitude',              'SemiMajor90', 'SemiMinor90', 'ErrorStrike', 'depth', 'depthError', 'magnitude', 'sigmaMagnitude','rake','dip','strike','type', 'probability',         'fuzzy']
        event=pandas.DataFrame(index=index,columns=columns)
        #event=pandas.Series()
        #assign event attributes
        event['eventID'] = smevent.attrib['event_id']
        event['Agency']  = smevent.attrib['event_network']
        event['year'],event['month'],event['day'],event['hour'],event['minute'],event['second'] = quakeml.utc2event(smevent.attrib['event_timestamp'])
        event['depth'] = float(smevent.attrib['depth'])
        event['magnitude'] = float(smevent.attrib['magnitude'])
        event['longitude'] = float(smevent.attrib['lon'])
        event['latitude'] = float(smevent.attrib['lat'])
        event['strike'] = float(smevent.attrib['strike'])
        event['dip'] = float(smevent.attrib['dip'])
        event['rake'] = float(smevent.attrib['rake'])
        event['type'] = shakeml.attrib['shakemap_event_type']
        #FIXME: deal with shakeml.attrib: 'shakemap_id': 'us1000gez7', 'shakemap_version': '2', 'code_version': '3.5.1615', 'process_timestamp': '2018-08-21T23:32:17Z', 'shakemap_originator': 'us', 'map_status': 'RELEASED'
        #FIXME: deal with description
        #smevent.attrib['event_description']

        #FIXME:uncertainty
        elems_event_specific_uncertainties = shakeml.findall('event_specific_uncertainty',namespaces = nsmap)
        index = [i for i in range(len(elems_event_specific_uncertainties))]
        columns = ['name','value','numsta']
        event_specific_uncertainties = pandas.DataFrame(index=index,columns=columns)
        for i,el in enumerate(elems_event_specific_uncertainties):
            event_specific_uncertainties.iloc[i]["name"] = el.attrib['name']
            event_specific_uncertainties.iloc[i].value = el.attrib['value']
            event_specific_uncertainties.iloc[i].numsta = el.attrib['numsta']

    #grid specification
    #NOTE:added indicator for structured and unstructured
    #TODO: derive regularity maybe...
    grid_specification=shakeml.find('grid_specification',namespaces = nsmap)
    try:
        regular_grid = bool(grid_specification.attrib['regular_grid'])
    except:
        #assume a regular grid
        regular_grid = True

    #TODO: actually necessary? Probably not...as is inherent to the grid if needed can be easily derived from pandas df
    #attributes: lon_min,lat_min,lon_max,lat_max,nominal_lon_spacing,nominal_lat_spacing,nlon,nlat

    #columns
    grid_fields = shakeml.findall('grid_field',namespaces = nsmap)

    #indices (start at 1) & argsort them
    column_idxs = [int(grid_field.attrib['index'])-1 for grid_field in grid_fields]
    idxs_sorted = np.argsort(column_idxs)
    column_names = [grid_field.attrib['name'] for grid_field in grid_fields]
    columns = [column_names[idx] for idx in idxs_sorted]

    #get grid
    grid_data = io.StringIO(shakeml.find('grid_data',namespaces = nsmap).text)

    grid_data = pandas.read_csv(grid_data,sep=' ',header=None)
    grid_data.columns = columns

    #get units
    units = pandas.DataFrame(index=[0],columns=columns)
    for grid_field in grid_fields:
        units.iloc[0][grid_field.attrib['name']] = grid_field.attrib['units']

    if siteml:
        return site_params(grid_data,units.iloc[0],regular_grid)
    else:
        #return a quakemap object
        return quakemap(event.iloc[0],event_specific_uncertainties,grid_data,units.iloc[0],regular_grid)

def quakemap2shakeml(qm,provider='GFZ'):
    '''
    Given a quakemap object generates a shakemap xml file (referred to as shakeml)
    Can also deal with a sites object
    '''
    try:
        event = qm.event
        event_specific_uncertainty = qm.event_specific_uncertainty
        shakemap = qm.shakemap
        siteml=False
    except:
        #Not a shakemap but a site map
        shakemap = qm.sites
        siteml=True

    units = qm.units
    regular_grid = qm.regular_grid

    #ensure that event is series
    #if type(event) != pandas.core.series.Series:
    #    print('WARNING: only implemented for one event, using first event of:\n{}'.format(event))
    #    event = event.iloc[0]

    nsmap = {'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
              None: 'http://earthquake.usgs.gov/eqcenter/shakemap'}
    schemaLocation = le.QName("{" + nsmap['xsi'] + "}schemaLocation")

    #processing attributes
    code_version = le.QName("code_version")
    shakemap_version = le.QName("shakemap_version")
    process_timestamp = le.QName("process_timestamp")
    shakemap_originator = le.QName("shakemap_originator")
    now = datetime.datetime.utcnow()
    now = pandas.Series({'year':now.year,'month':now.month,'day':now.day,'hour':now.hour,'minute':now.minute,'second':now.second+now.microsecond/10.**6})
    if siteml:
        shakeml = le.Element('site_grid',
                            {schemaLocation:"http://earthquake.usgs.gov http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd",
                             #event_id: event.eventID,
                             #FIXME: same as eventID!? No should be related to measure, gmpe etc....
                             #shakemap_id: event.eventID,
                             #NOTE: not shakemap standard
                             code_version: 'shakyground 0.1',
                             shakemap_version: '1',
                             process_timestamp: quakeml.event2utc(now),
                             shakemap_originator: provider,
                             #map_status: 'RELEASED',
                             #shakemap_event_type: event.type,
                            }, nsmap=nsmap)
    else:
        event_id    = le.QName("event_id")
        shakemap_id = le.QName("shakemap_id")
        map_status = le.QName("map_status")
        shakemap_event_type = le.QName("shakemap_event_type")
        shakeml = le.Element('shakemap_grid',
                            {schemaLocation:"http://earthquake.usgs.gov http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd",
                             event_id: event.eventID,
                             #FIXME: same as eventID!? No should be related to measure, gmpe etc....
                             shakemap_id: event.eventID,
                             #NOTE: not shakemap standard
                             code_version: 'shakyground 0.1',
                             shakemap_version: '1',
                             process_timestamp: quakeml.event2utc(now),
                             shakemap_originator: provider,
                             map_status: 'RELEASED',
                             shakemap_event_type: event.type,
                            }, nsmap=nsmap)

        # write event data
        # <event event_id="us1000gez7" magnitude="7.3" depth="123.18" lat="10.739200" lon="-62.910600" event_timestamp="2018-08-21T21:31:42UTC"                  event_network="us" event_description="OFFSHORE SUCRE, VENEZUELA" />
        magnitude = le.QName("magnitude")
        depth = le.QName("depth")
        lat = le.QName("lat")
        lon = le.QName("lon")
        strike = le.QName("strike")
        rake = le.QName("rake")
        dip = le.QName("dip")
        event_timestamp = le.QName("event_timestamp")
        event_network = le.QName("event_network")
        event_description = le.QName("event_description")
        smevent = le.SubElement(shakeml,'event',
                               {event_id: str(event.eventID),
                                magnitude: str(event.magnitude),
                                depth: str(event.depth),
                                lat:str(event.latitude),
                                lon:str(event.longitude),
                                strike:str(event.strike),
                                rake:str(event.rake),
                                dip:str(event.dip),
                                event_timestamp:str(quakeml.event2utc(event)),
                                event_network:str(event.Agency),
                                event_description:''
                                },
                                nsmap=nsmap
                               )

    # write metadata on grid
    # <grid_specification lon_min="-67.910600" lat_min="5.829200" lon_max="-57.910600" lat_max="15.649200" nominal_lon_spacing="0.016667" nominal_lat_spacing="0.016672" nlon="601" nlat="590" />
    lon_min = le.QName("lon_min")
    lat_min = le.QName("lat_min")
    lon_max = le.QName("lon_max")
    lat_max = le.QName("lat_max")
    nominal_lon_spacing = le.QName("nominal_lon_spacing")
    nominal_lat_spacing = le.QName("nominal_lat_spacing")
    nlon = le.QName("nlon")
    nlat = le.QName("nlat")
    reg_grid = le.QName("regular_grid")
    # get plon and plat
    if regular_grid:
        grid_specification = le.SubElement(shakeml, 'grid_specification',
                                           {lon_min: str(shakemap.LON.min()),
                                            lat_min: str(shakemap.LAT.min()),
                                            lon_max: str(shakemap.LON.max()),
                                            lat_max: str(shakemap.LAT.max()),
                                            nominal_lon_spacing: str(
                                                round(abs(np.mean(np.diff(shakemap.LON.unique())[:-1])), 6)),
                                            nominal_lat_spacing: str(
                                                round(abs(np.mean(np.diff(shakemap.LAT.unique())[:-1])), 6)),
                                            nlon: str(len(shakemap.LON.unique())),
                                            nlat: str(len(shakemap.LAT.unique())),
                                            reg_grid: '1'
                                            },
                                           nsmap=nsmap
                                           )
    else:
        grid_specification = le.SubElement(shakeml, 'grid_specification',
                                           {lon_min: str(shakemap.LON.min()),
                                            lat_min: str(shakemap.LAT.min()),
                                            lon_max: str(shakemap.LON.max()),
                                            lat_max: str(shakemap.LAT.max()),
                                            reg_grid: '0'
                                            },
                                           nsmap=nsmap
                                           )

    if not siteml:
        #FIXME: which use somewhere on our side??
        #<event_specific_uncertainty name="pga" value="0.000000" numsta="" />
        #<event_specific_uncertainty name="pgv" value="0.000000" numsta="" />
        #<event_specific_uncertainty name="mi" value="0.000000" numsta="" />
        #<event_specific_uncertainty name="psa03" value="0.000000" numsta="" />
        #<event_specific_uncertainty name="psa10" value="0.000000" numsta="" />
        #<event_specific_uncertainty name="psa30" value="0.000000" numsta="" />
        list_event_specific_uncertainty=[]
        name = le.QName("name")
        value = le.QName("value")
        numsta = le.QName("numsta")
        for i in range(len(event_specific_uncertainty)):#["pga","pgv","mi","psa03","psa10","psa30"]:
            list_event_specific_uncertainty.append(le.SubElement(shakeml,'event_specific_uncertainty',
                                                      {name: str(event_specific_uncertainty.iloc[i]["name"]),
                                                       value: str(event_specific_uncertainty.iloc[i]["value"]),
                                                       numsta:str(event_specific_uncertainty.iloc[i]["numsta"])},
                                                      nsmap=nsmap))

    #grid field specification
    #<grid_field index="1" name="LON" units="dd" />
    #<grid_field index="2" name="LAT" units="dd" />
    #<grid_field index="3" name="PGA" units="pctg" />
    #<grid_field index="4" name="PGV" units="cms" />
    #<grid_field index="5" name="MMI" units="intensity" />
    #<grid_field index="6" name="PSA03" units="pctg" />
    #<grid_field index="7" name="PSA10" units="pctg" />
    #<grid_field index="8" name="PSA30" units="pctg" />
    #<grid_field index="9" name="STDPGA" units="ln(pctg)" />
    #<grid_field index="10" name="URAT" units="" />
    #<grid_field index="11" name="SVEL" units="ms" />
    index = le.QName("index")
    _name = le.QName("name")
    _units = le.QName("units")
    grid_fields=[]
    for i,col in enumerate(shakemap.columns):
        grid_fields.append(le.SubElement(shakeml,'grid_field',
                                        {index:str(i+1),#starts at 1
                                         _name: col,
                                         _units: str(units[col])},
                                         nsmap=nsmap))

    #grid data
    grid_data=le.SubElement(shakeml,'grid_data',nsmap=nsmap)
    grid_data.text = '\n'+shakemap.to_csv(sep=' ',header=False, index=False)
    #grid_data.text = '\n'+shakemap.to_string(header=False,index=False,justify='left')

    return le.tostring(shakeml,pretty_print=True,encoding='unicode')
