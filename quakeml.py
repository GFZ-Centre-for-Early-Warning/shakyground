#####################################
# Convert quakeml catalogs to pandas
# and vice versa
import pandas
import lxml.etree as le

#TODO: publicID in quakeml should refer to webservice address

def event2utc(event):
    '''
    given event returns UTC string
    '''
    d=event.fillna(0)
    return '{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:09f}Z'.format(int(d.year),int(max(d.month,1)),int(max(d.day,1)),int(d.hour),int(d.minute),d.second)

def utc2event(utc):
    '''
    given utc string returns list with year,month,day,hour,minute,second
    '''
    #last part usually either Z(ulu) or UTC, if not fails
    if utc[-3:]=='UTC':
        utc=utc[:-2]
    elif utc[-1:]=='Z':
        pass
    else:
        raise Exception('Cannot handle timezone other than Z(ulu) or UTC: {}'.format(utc))
    date,time = utc.split('T')
    return [int(v) if i<5 else float(v) for i,v in enumerate([int(d) for d in date.split('-')]+[float(t) for t in time[:-1].split(':')])]

def events2quakeml(catalog,provider='GFZ'):
    '''
    Given a pandas dataframe with events returns QuakeML version of
    the catalog
    '''
    #TODO: add uncertainty to all values (NOTE: OQ/HMTK style spatial uncertainty is ellipse semi-major/semi-minor/strike-error)
    xml_namespace = 'http://quakeml.org/xmlns/quakeml/1.2'
    quakeml = le.Element('eventParameters',namespace=xml_namespace)
    #go through all events
    for i in range(len(catalog)):
        quake = catalog.iloc[i]
        event = le.SubElement(quakeml,'event',{'publicID':str(quake.eventID)})
        preferredOriginID = le.SubElement(event,'preferredOriginID')
        preferredOriginID.text=str(quake.eventID)
        preferredMagnitudeID = le.SubElement(event,'preferredMagnitudeID')
        preferredMagnitudeID.text=str(quake.eventID)
        qtype = le.SubElement(event,'type')
        qtype.text = 'earthquake'
        description = le.SubElement(event,'description')
        text = le.SubElement(description,'text')
        text.text = str(quake.type)
        #origin
        origin = le.SubElement(event,'origin',{'publicID':str(quake.eventID)})
        time = le.SubElement(origin,'time')
        value = le.SubElement(time,'value')
        value.text = event2utc(quake)
        latitude = le.SubElement(origin,'latitude')
        value = le.SubElement(latitude,'value')
        value.text = str(quake.latitude)
        longitude = le.SubElement(origin,'longitude')
        value = le.SubElement(longitude,'value')
        value.text = str(quake.longitude)
        depth = le.SubElement(origin,'depth')
        value = le.SubElement(depth,'value')
        value.text = str(quake.depth)
        creationInfo = le.SubElement(origin,'creationInfo')
        author = le.SubElement(creationInfo,'value')
        author.text = provider
        #magnitude
        magnitude = le.SubElement(event,'magnitude',{'publicID':str(quake.eventID)})
        mag = le.SubElement(magnitude,'mag')
        value = le.SubElement(mag,'value')
        value.text = str(quake.magnitude)
        mtype = le.SubElement(magnitude,'type')
        mtype.text = 'MW'
        creationInfo = le.SubElement(magnitude,'creationInfo')
        author = le.SubElement(creationInfo,'value')
        author.text = provider
        #plane (write only fault plane not auxilliary)
        focalMechanism = le.SubElement(event,'focalMechanism',{'publicID':str(quake.eventID)})
        nodalPlanes = le.SubElement(focalMechanism,'nodalPlanes')
        nodalPlane1 = le.SubElement(nodalPlanes,'nodalPlane1')
        strike = le.SubElement(nodalPlane1,'strike')
        value  = le.SubElement(strike,'value')
        value.text = str(quake.strike)
        dip = le.SubElement(nodalPlane1,'dip')
        value  = le.SubElement(dip,'value')
        value.text = str(quake.dip)
        rake = le.SubElement(nodalPlane1,'rake')
        value  = le.SubElement(rake,'value')
        value.text = str(quake.rake)
        preferredPlane = le.SubElement(nodalPlanes,'preferredPlane')
        preferredPlane.text = 'nodalPlane1'

    #return str(le.tostring(quakeml,pretty_print=True,xml_declaration=True),encoding='utf-8')
    return le.tostring(quakeml,pretty_print=True,encoding='unicode')

def quakeml2events(quakemlfile,provider='GFZ'):
    '''
    Given a quakeml file/or string returns a pandas dataframe
    '''
    #TODO: add uncertainty
    try:
        #read quakeml catalog
        with open(quakemlfile,'r') as f:
            quakeml = f.read()
    except:
        #maybe already string
        quakeml = quakemlfile

    quakeml = le.fromstring(quakeml)
    #initialize catalog
    index = [i for i in range(len(quakeml))]
    columns=['eventID', 'Agency', 'Identifier', 'year', 'month', 'day', 'hour', 'minute', 'second', 'timeError', 'longitude', 'latitude',              'SemiMajor90', 'SemiMinor90', 'ErrorStrike', 'depth', 'depthError', 'magnitude', 'sigmaMagnitude','rake','dip','strike','type', 'probability',   'fuzzy']
    #columns=['eventID', 'Agency', 'Identifier', 'year', 'month', 'day', 'hour', 'minute', 'second', 'timeError', 'longitude', 'latitude','SemiMajor90', 'SemiMinor90', 'ErrorStrike', 'depth', 'depthError', 'magnitude', 'sigmaMagnitude', 'type', 'probability', 'fuzzy']
    catalog=pandas.DataFrame(index=index,columns=columns)
    #add individual events to catalog
    for i,event in enumerate(quakeml):
        #get ID
        catalog.iloc[i].eventID = event.attrib['publicID']
        #type
        catalog.iloc[i].type = event.find('description').findtext('text')
        #origin
        origin = event.find('origin')
        #time
        catalog.iloc[i].year,catalog.iloc[i].month,catalog.iloc[i].day,catalog.iloc[i].hour,catalog.iloc[i].minute,catalog.iloc[i].second = utc2event(origin.find('time').findtext('value'))
        #latitude/longitude/depth
        catalog.iloc[i].latitude = float(origin.find('latitude').findtext('value'))
        catalog.iloc[i].longitude = float(origin.find('longitude').findtext('value'))
        catalog.iloc[i].depth = float(origin.find('depth').findtext('value'))
        #agency/provider
        catalog.iloc[i].agency = origin.find('creationInfo').findtext('value')
        #magnitude
        catalog.iloc[i].magnitude = float(event.find('magnitude').find('mag').findtext('value'))
        #plane
        nodalPlanes = event.find('focalMechanism').find('nodalPlanes')
        preferredPlane = nodalPlanes.findtext('preferredPlane')
        catalog.iloc[i].strike = float(nodalPlanes.find(preferredPlane).find('strike').findtext('value'))
        catalog.iloc[i].dip = float(nodalPlanes.find(preferredPlane).find('dip').findtext('value'))
        catalog.iloc[i].rake = float(nodalPlanes.find(preferredPlane).find('rake').findtext('value'))

    return catalog
