==========================
Boston Scientific Cartesia
==========================

* **Manufacturer:** `Boston Scientific <https://www.bostonscientific.com/en-US/home.html>`_ 
* **Product:** Boston Scientific Cartesia  
* **Serial name:** CartesiaX, CartesiaHX  

Source documentation:  
`Source PDF <https://www.bostonscientific.com/content/dam/elabeling/nm/92495783-02_Vercis_TM-DBS_Systems_Surgical_Implant_Manual_multi-OUS_s.pdf>`_  
`Source Webpage <https://www.bostonscientific.com/en-US/products/deep-brain-stimulation-systems/vercise-genus-dbs-system/leads.html>`_  

The source documentation also contains information about :doc:`Boston Scientific Directed Lead (DB-2202) <./Boston_Scientific_Vercise_Directed>`.

---------------------------
Boston Scientific CartesiaX
---------------------------

~~~~~~~~~~~~~~~~~~~~~~~
Default Parameters (mm)
~~~~~~~~~~~~~~~~~~~~~~~

* tip_length = 1.1
* contact_length = 1.5
* contact_spacing = 0.5
* lead_diameter = 1.3
* total_length = 450.0
* contacts_skipped = 5.0

----
Code
----

.. autoclass:: ossdbs.electrodes.boston_scientific_cartesia.BostonScientificCartesiaXModel
    :members:  
    :show-inheritance:  

----------------------------
Boston Scientific CartesiaHX
----------------------------

~~~~~~~~~~~~~~~~~~~~~~~
Default Parameters (mm)
~~~~~~~~~~~~~~~~~~~~~~~

* tip_length = 1.1
* contact_length = 1.5
* contact_spacing = 0.5
* lead_diameter = 1.3
* total_length = 450.0
* contacts_skipped = 7.0

----
Code
----

.. autoclass:: ossdbs.electrodes.boston_scientific_cartesia.BostonScientificCartesiaHXModel
    :members:  
    :show-inheritance:  
