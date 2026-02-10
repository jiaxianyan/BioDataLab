import requests
import json

def get_xml_by_ndc(ndc_code):
    base_url = "https://dailymed.nlm.nih.gov/dailymed/services/v2"
    
    search_url = f"{base_url}/spls.json"
    params = {'ndc': ndc_code}
    
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        
        data = response.json()
        
        if not data.get('data'):
            return None
            
        first_result = data['data'][0]
        set_id = first_result.get('setid')
        
    except Exception as e:
        return None

    xml_url = f"{base_url}/spls/{set_id}.xml"
    
    try:
        xml_response = requests.get(xml_url)
        xml_response.raise_for_status()
        
        return xml_response.text
        
    except Exception as e:
        return None

if __name__ == "__main__":
    target_ndc = "42681-0025-1"
    
    xml_content = get_xml_by_ndc(target_ndc)
    
    if xml_content:
        with open(f"{target_ndc}.xml", "w", encoding="utf-8") as f:
            f.write(xml_content)