

 ////copy paste this scrip in the same directory as PhamCat Class

 package org.pharmgkb.pharmcat;

 import java.io.File;
 import java.io.FilenameFilter;
 
 
 public class Runerz
 {
   public static void main(String[] args) {
 //    String directory = "/home/xixe/pharm/pharmvip-guideline/resources/samples/bigchunk";
     String directory = "/home/xixe/pharm/pharmvip-guideline/resources/samples/mini_chunk/390";
 //    File folder = new File("/home/xixe/pharm/pharmvip-guideline/resources/samples/mini_chunk/");
 //    File[] listOfFiles = folder.listFiles();
 //
 //    for (int i = 0; i < listOfFiles.length; i++) {
 //      if (listOfFiles[i].isFile()) {
 //        System.out.println("File " + listOfFiles[i].getName());
 //      } else if (listOfFiles[i].isDirectory()) {
 //        System.out.println("Directory " + listOfFiles[i].getName());
 //      }
 //    }
     File dir = new File(directory);
 
     File[] matches = dir.listFiles(new FilenameFilter()
     {
       public boolean accept(File dir, String name)
       {
 //        return name.startsWith("temp") && name.endsWith(".vcf");
         return name.endsWith(".vcf");
       }
     });
 
     for (int i = 0; i < matches.length; i++) {
       System.out.println(matches[i]);
     }
 
     for (int i = 0; i < matches.length; i++) {
       String[] arguments = new String[] {"-vcf",matches[i].toString(),
                 "-o","/home/xixe/tmp/optt/phCat"};
       PharmCAT.main(arguments);
     }
 
 
 //
   }
 }
 