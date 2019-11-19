import os
import jinja2

def get_templates():
    current_location = os.path.dirname(os.path.abspath(__file__))
    for root, dirs, files in os.walk(current_location):
        if "templates" in dirs:
            return os.path.join(root, "templates")

templateLoader = jinja2.FileSystemLoader(searchpath=get_templates())
templateEnv = jinja2.Environment(loader=templateLoader)

def mk_gal(images, outfile):

    TEMPLATE_FILE = "gallery_lazy.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'ncol':4, 'width': 230, 'height': 225}
    outputText = template.render(imlist=images, opt=opt)  # this is where to put args to the template renderer

    #print(outputText)

    with open(outfile, "wb") as fh:
        fh.write(outputText)

def mk_cohort_gal(section_list, outfile):

    TEMPLATE_FILE = "cohort_viz_gallery_lazy.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'ncol':4, 'width': 230, 'height': 225}

    outputText = template.render(sclist=section_list, opt=opt)  # this is where to put args to the template renderer

    #print(outputText)

    with open(outfile, "wb") as fh:
        fh.write(outputText)

def mk_index(table_headers, table_tuples, outfolder, project_name):

    TEMPLATE_FILE = "index.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    opt = {'title':project_name}
    outputText = template.render(headers=table_headers, summary=table_tuples, opt=opt)

    with open(os.path.join(outfolder, "index.html"), "wb") as fh:
        fh.write(outputText)
    return True

