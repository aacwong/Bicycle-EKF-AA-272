////////////////////////////////////////////////////////////////////////////////
// The following FIT Protocol software provided may be used with FIT protocol
// devices only and remains the copyrighted property of Garmin Canada Inc.
// The software is being provided on an "as-is" basis and as an accommodation,
// and therefore all warranties, representations, or guarantees of any kind
// (whether express, implied or statutory) including, without limitation,
// warranties of merchantability, non-infringement, or fitness for a particular
// purpose, are specifically disclaimed.
//
// Copyright 2021 Garmin International, Inc.
////////////////////////////////////////////////////////////////////////////////
// ****WARNING****  This file is auto-generated!  Do NOT edit this file.
// Profile Version = 21.67Release
// Tag = production/akw/21.67.00-0-gd790f76b
////////////////////////////////////////////////////////////////////////////////


package com.garmin.fit.csv;

import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import com.garmin.fit.*;

public class CSVWriter {
    private String fileName;
    private File file;
    private BufferedWriter writer;
    private ArrayList<String> headers;
    private ArrayList<String> values;

    public CSVWriter(String fileName) {
        this.fileName = fileName;
        this.headers = new ArrayList<String>();
        this.values = new ArrayList<String>();
    }

    public void close() {
        try {
            if (writer != null) {
                BufferedReader reader;

                writer.close();
                writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), "UTF-8"));
                // Adds the UTF-8 BOM so editors can recognize the file as UTF-8
                writer.write(Fit.UTF8_BOM);

                for (int i = 0; i < headers.size(); i++) {
                    writer.write(headers.get(i) + ",");
                }

                writer.write("\n");

                reader = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"));

                while (reader.ready()) {
                    writer.write(reader.readLine() + "\n");
                }

                reader.close();
                writer.close();
                file.delete();
            }
        } catch (java.io.IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void clear() {
        for (int i = 0; i < values.size(); i++) {
            values.set(i, new String(""));
        }
    }

    public void set(String header, Object value) {
        if (header == null) {
            header = "null";
        }

        if (value == null) {
            value = "null";
        }

        for (int i = 0; i < headers.size(); i++) {
            if (headers.get(i).compareTo(header) == 0) {
                values.set(i, value.toString());
                return;
            }
        }

        headers.add(header.toString());
        values.add(value.toString());
    }

    public void writeln() {
        try {
            if (writer == null) {
                file = new File(fileName + ".tmp");
                writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "UTF-8"));
            }

            for (int i = 0; i < values.size(); i++) {
                writer.write(values.get(i) + ",");
            }

            writer.write("\n");
        } catch (java.io.IOException e) {
            throw new RuntimeException(e);
        }
    }
}
