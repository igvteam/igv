package org.igv.feature.genome;

import org.json.JSONObject;

public class ClinVar {

    public static String getClinVarURL(String hgvsNotation) {
        try {
            String encodedHgvs = java.net.URLEncoder.encode(hgvsNotation, java.nio.charset.StandardCharsets.UTF_8);
            String esearchUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" +
                    "db=clinvar&term=" + encodedHgvs + "&retmode=json";

            java.net.http.HttpClient client = java.net.http.HttpClient.newHttpClient();
            java.net.http.HttpRequest request = java.net.http.HttpRequest.newBuilder()
                    .uri(java.net.URI.create(esearchUrl))
                    .build();

            java.net.http.HttpResponse<String> response = client.send(request,
                    java.net.http.HttpResponse.BodyHandlers.ofString());

            // Parse JSON response to get the first ClinVar accession
            String body = response.body();

            JSONObject json = new JSONObject(body);
            JSONObject esearchResult = json.getJSONObject("esearchresult");
            if (esearchResult.getInt("count") > 0) {
                String uid = esearchResult.getJSONArray("idlist").getString(0);
                return "https://www.ncbi.nlm.nih.gov/clinvar/variation/" + uid  + "/";
            } else {
                return null;
            }

        } catch (Exception e) {
            return null;
        }
    }
}
