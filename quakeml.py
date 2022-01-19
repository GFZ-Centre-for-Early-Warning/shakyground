#####################################
# Convert quakeml catalogs to pandas
# and vice versa
import pandas
import lxml.etree as le

# example value, will be set later by quakeml2events
NS_QUAKEML = "http://quakeml.org/xmlns/bed/1.2"

# TODO: publicID in quakeml should refer to webservice address


def event2utc(event):
    """
    given event returns UTC string
    """
    d = event.fillna(0)
    return "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:09f}Z".format(
        int(d.year),
        int(max(d.month, 1)),
        int(max(d.day, 1)),
        int(d.hour),
        int(d.minute),
        d.second,
    )


def utc2event(utc):
    """
    given utc string returns list with year,month,day,hour,minute,second
    """
    # last part usually either Z(ulu) or UTC, if not fails
    if utc[-3:] == "UTC":
        utc = utc[:-2]
    elif utc[-1:] == "Z":
        pass
    else:
        raise Exception(
            "Cannot handle timezone other than Z(ulu) or UTC: {}".format(utc)
        )
    date, time = utc.split("T")
    return [
        int(v) if i < 5 else float(v)
        for i, v in enumerate(
            [int(d) for d in date.split("-")]
            + [float(t) for t in time[:-1].split(":")]
        )
    ]


def events2quakeml(catalog, provider="GFZ"):
    """
    Given a pandas dataframe with events returns QuakeML version of
    the catalog
    """
    # TODO: add uncertainty to all values (NOTE: OQ/HMTK style spatial uncertainty is ellipse semi-major/semi-minor/strike-error)
    quakeml = le.Element("eventParameters", namespace=NS_QUAKEML)
    # go through all events
    for i in range(len(catalog)):
        quake = catalog.iloc[i]
        event = le.SubElement(
            quakeml, "event", {"publicID": str(quake.eventID)}
        )
        preferredOriginID = le.SubElement(
            event, "preferredOriginID", namespace=NS_QUAKEML
        )
        preferredOriginID.text = str(quake.eventID)
        preferredMagnitudeID = le.SubElement(
            event, "preferredMagnitudeID", namespace=NS_QUAKEML
        )
        preferredMagnitudeID.text = str(quake.eventID)
        qtype = le.SubElement(event, "type", namespace=NS_QUAKEML)
        qtype.text = "earthquake"
        description = le.SubElement(event, "description", namespace=NS_QUAKEML)
        text = le.SubElement(description, "text", namespace=NS_QUAKEML)
        text.text = str(quake.type)
        # origin
        origin = le.SubElement(
            event,
            "origin",
            {"publicID": str(quake.eventID)},
            namespace=NS_QUAKEML,
        )
        time = le.SubElement(origin, "time", namespace=NS_QUAKEML)
        value = le.SubElement(time, "value", namespace=NS_QUAKEML)
        value.text = event2utc(quake)
        latitude = le.SubElement(origin, "latitude", namespace=NS_QUAKEML)
        value = le.SubElement(latitude, "value", namespace=NS_QUAKEML)
        value.text = str(quake.latitude)
        longitude = le.SubElement(origin, "longitude", namespace=NS_QUAKEML)
        value = le.SubElement(longitude, "value", namespace=NS_QUAKEML)
        value.text = str(quake.longitude)
        depth = le.SubElement(origin, "depth", namespace=NS_QUAKEML)
        value = le.SubElement(depth, "value", namespace=NS_QUAKEML)
        value.text = str(quake.depth)
        creationInfo = le.SubElement(
            origin, "creationInfo", namespace=NS_QUAKEML
        )
        author = le.SubElement(creationInfo, "value", namespace=NS_QUAKEML)
        author.text = provider
        # magnitude
        magnitude = le.SubElement(
            event,
            "magnitude",
            {"publicID": str(quake.eventID)},
            namespace=NS_QUAKEML,
        )
        mag = le.SubElement(magnitude, "mag", namespace=NS_QUAKEML)
        value = le.SubElement(mag, "value", namespace=NS_QUAKEML)
        value.text = str(quake.magnitude)
        mtype = le.SubElement(magnitude, "type", namespace=NS_QUAKEML)
        mtype.text = "MW"
        creationInfo = le.SubElement(
            magnitude, "creationInfo", namespace=NS_QUAKEML
        )
        author = le.SubElement(creationInfo, "value", namespace=NS_QUAKEML)
        author.text = provider
        # plane (write only fault plane not auxilliary)
        focalMechanism = le.SubElement(
            event,
            "focalMechanism",
            {"publicID": str(quake.eventID)},
            namespace=NS_QUAKEML,
        )
        nodalPlanes = le.SubElement(
            focalMechanism, "nodalPlanes", namespace=NS_QUAKEML
        )
        nodalPlane1 = le.SubElement(
            nodalPlanes, "nodalPlane1", namespace=NS_QUAKEML
        )
        strike = le.SubElement(nodalPlane1, "strike", namespace=NS_QUAKEML)
        value = le.SubElement(strike, "value", namespace=NS_QUAKEML)
        value.text = str(quake.strike)
        dip = le.SubElement(nodalPlane1, "dip", namespace=NS_QUAKEML)
        value = le.SubElement(dip, "value", namespace=NS_QUAKEML)
        value.text = str(quake.dip)
        rake = le.SubElement(nodalPlane1, "rake", namespace=NS_QUAKEML)
        value = le.SubElement(rake, "value", namespace=NS_QUAKEML)
        value.text = str(quake.rake)
        preferredPlane = le.SubElement(
            nodalPlanes, "preferredPlane", namespace=NS_QUAKEML
        )
        preferredPlane.text = "nodalPlane1"

    # return str(le.tostring(quakeml,pretty_print=True,xml_declaration=True),encoding='utf-8')
    return le.tostring(quakeml, pretty_print=True, encoding="unicode")


def find_element(element, which):
    """
    :type element: xml.etree.ElementTree.Element
    :type which: str
    :rtype: xml.etree.ElementTree.Element
    """
    return element.find("{" + NS_QUAKEML + "}" + which)


def find_text(element, text):
    """
    :type element: xml.etree.ElementTree.Element
    :type text: str
    :rtype: str
    """
    return element.findtext("{" + NS_QUAKEML + "}" + text)


def find_text_in_element(element, which, text):
    elm = find_element(element, which)

    if elm is None:
        return None

    return find_text(elm, text)


def quakeml2events(quakemlfile, provider="GFZ"):
    """
    Given a quakeml file/or string returns a pandas dataframe
    """
    # TODO: add uncertainty
    try:
        # read quakeml catalog
        with open(quakemlfile, "r") as f:
            quakeml = f.read()
    except:
        # maybe already string
        quakeml = quakemlfile

    quakeml = le.fromstring(quakeml)

    # initialize catalog
    index = [i for i in range(len(quakeml))]
    columns = [
        "eventID",
        "Agency",
        "Identifier",
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "timeError",
        "longitude",
        "latitude",
        "SemiMajor90",
        "SemiMinor90",
        "ErrorStrike",
        "depth",
        "depthError",
        "magnitude",
        "sigmaMagnitude",
        "rake",
        "dip",
        "strike",
        "type",
        "probability",
        "fuzzy",
    ]
    # columns=['eventID', 'Agency', 'Identifier', 'year', 'month', 'day', 'hour', 'minute', 'second', 'timeError', 'longitude', 'latitude','SemiMajor90', 'SemiMinor90', 'ErrorStrike', 'depth', 'depthError', 'magnitude', 'sigmaMagnitude', 'type', 'probability', 'fuzzy']
    catalog = pandas.DataFrame(index=index, columns=columns)

    # add individual events to catalog
    for i, event in enumerate(quakeml):
        # get ID
        catalog.iloc[i].eventID = event.attrib["publicID"]
        # type
        catalog.iloc[i].type = find_text_in_element(
            event, "description", "text"
        )
        # origin
        origin = find_element(event, "origin")
        # time
        timevalue = find_text_in_element(origin, "time", "value")
        (
            catalog.iloc[i].year,
            catalog.iloc[i].month,
            catalog.iloc[i].day,
            catalog.iloc[i].hour,
            catalog.iloc[i].minute,
            catalog.iloc[i].second,
        ) = utc2event(timevalue)
        # latitude/longitude/depth
        catalog.iloc[i].latitude = float(
            find_text_in_element(origin, "latitude", "value")
        )
        catalog.iloc[i].longitude = float(
            find_text_in_element(origin, "longitude", "value")
        )
        catalog.iloc[i].depth = float(
            find_text_in_element(origin, "depth", "value")
        )
        # agency/provider
        catalog.iloc[i].agency = find_text_in_element(
            origin, "creationInfo", "value"
        )
        # magnitude
        mag = find_element(event, "magnitude")
        catalog.iloc[i].magnitude = float(
            find_text_in_element(mag, "mag", "value")
        )
        # plane
        focal = find_element(event, "focalMechanism")
        nodalPlanes = find_element(focal, "nodalPlanes")
        preferredPlane = nodalPlanes.attrib["preferredPlane"]

        if preferredPlane == "1":
            preferredPlane = "nodalPlane1"
        elif preferredPlane == "1":
            preferredPlane = "nodalPlane2"

        preferredPlaneElm = find_element(nodalPlanes, preferredPlane)

        if preferredPlaneElm is None:
            continue

        catalog.iloc[i].strike = float(
            find_text_in_element(preferredPlaneElm, "strike", "value")
        )
        catalog.iloc[i].dip = float(
            find_text_in_element(preferredPlaneElm, "dip", "value")
        )
        catalog.iloc[i].rake = float(
            find_text_in_element(preferredPlaneElm, "rake", "value")
        )

    return catalog
