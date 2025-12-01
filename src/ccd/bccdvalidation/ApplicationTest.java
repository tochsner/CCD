package ccd.bccdvalidation;

import beast.base.inference.Runnable;
import beastfx.app.tools.Application;

public class ApplicationTest extends Runnable {
    public static void main(String[] args) throws Exception {
        new Application(new ApplicationTest(), "Application Test", args);
    }

    @Override
    public void initAndValidate() {

    }

    @Override
    public void run() throws Exception {
        System.out.println("Hi!");
    }
}
