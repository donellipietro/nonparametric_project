import http.client, urllib.request, urllib.parse, urllib.error, base64

headers = {
    # Request headers
    'Ocp-Apim-Subscription-Key': '404c559bf2bf48ce8ec111cd72bfcfe4',
}

params = urllib.parse.urlencode({
})

try:
    conn = http.client.HTTPSConnection('api.delijn.be')
    conn.request("GET", "/DLKernOpenData/v1/beta/haltes?%s" % params, "{body}", headers)
    response = conn.getresponse()
    data = response.read()
    conn.close()
except Exception as e:
    print("[Errno {0}] {1}".format(e.errno, e.strerror))

import json
JSON_NAME = 'data/output.json'
data_s = data.decode('utf-8')
json_data = json.loads(data_s)
with open(JSON_NAME, 'w') as another_open_file:
    another_open_file.write(data_s)