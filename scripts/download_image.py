import mechanize
from BeautifulSoup import BeautifulSoup as BS
import matplotlib.pyplot as plt
import pyfits as pf
from astLib import astWCS, astPlots

RA = '05 22 59.0' 
Dec = '-36 27 30.2'
# Browser
br = mechanize.Browser()


# Browser options - need these for good browser behaviour?
br.set_handle_equiv(True)
#br.set_handle_gzip(True)
br.set_handle_redirect(True)
br.set_handle_referer(True)
br.set_handle_robots(False)

r = br.open('http://www.cv.nrao.edu/nvss/postage.shtml')

#print br.title()


##Find out the name of all forms on the website
#for form in br.forms():
	#print form.name
	#print form
	
##Form names print as None type, so use this to select the first form:
br.select_form(nr=0)

#for control in br.form.controls:
    #print control
    #print "type=%s, name=%s value=%s" % (control.type, control.name, br[control.name])

br.form['RA'] = RA
br.form['Dec'] = Dec

br['Type'] = ['application/octet-stream']

#for form in br.forms(): print form

for control in br.form.controls:
   if control.type == "submit":
       control.disabled = True

response = br.submit()

#print type(response.read())

soup = BS(response.read())

body_tag = soup.body













#br.submit()



#print br.form

#response = br.response().read()

#print response

#with open("mechanize_results.html", "w") as f:
	#f.write(response)

#print response.read()



#print type(soup)

 
#browser.open('http://www.cv.nrao.edu/nvss/postage.shtml')
#browser.select_form(name='the_form')
#browser['field1'] = 'value'
#browser['field2'] = 'value'
#browser['field3'] = 'value'
#browser.submit()


