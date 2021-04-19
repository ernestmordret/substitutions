from lxml import etree
import pandas as pd


# Create a tree for the XML document
with open('./unimod.xml') as f:
    tree = etree.parse(f)


def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]

conflicting_classes = ['Post-translational','Artefact', 'Pre-translational',
                       'Other glycosylation', 'O-linked glycosylation',
                       'N-linked glycosylation',
                       'Non-standard residue', 'Co-translational']


root = tree.getroot()
for ns in root.nsmap.values():
    remove_namespace(tree, ns)

modifications = root.find('modifications')

lines = []

for mod in modifications:
    title = mod.attrib.get('title')
    base_delta = float(mod.find('delta').attrib.get('mono_mass'))
    for spec in mod.findall('specificity'):
        site = spec.attrib.get('site')
        position = spec.attrib.get('position')
        classification = spec.attrib.get('classification')
        
        if spec.find('NeutralLoss') is not None:
            for neutral_loss in spec.findall('NeutralLoss'):

                comp = neutral_loss.attrib.get('composition')
                neutral_loss_mass = float(neutral_loss.attrib.get('mono_mass'))
                lines.append({'title': f'{title} - {comp}',
                             'delta_m' : base_delta - neutral_loss_mass,
                             'classification' : classification,
                             'site': site,
                             'position': position})
        else:
            
            lines.append({'title': title,
                         'delta_m' : base_delta,
                         'classification' : classification,
                         'site': site,
                         'position': position})

mods = pd.DataFrame(lines)
subs = mods[mods['classification']=='AA substitution']
mods = mods[mods['classification'].isin(conflicting_classes) & (mods.delta_m!=0)]
mods = mods.sort_values('delta_m')
mods.reset_index(drop=True).to_pickle('./danger_mods')

danger_mods = pd.read_pickle('./danger_mods')