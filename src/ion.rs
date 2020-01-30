use ptable::Element;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Ion {
    element: Element,
    charge: i8
}

impl Ion {
    pub fn new(element: Element, charge: i8) -> Ion {
        Ion { element, charge }
    }

    #[inline(always)]
    pub fn get_element(&self) -> &Element {
        &self.element
    }

    #[inline(always)]
    pub fn get_element_mut(&mut self) -> &mut Element {
        &mut self.element
    }

    #[inline(always)]
    pub fn set_element(&mut self, element: Element) {
        self.element = element;
    }

    #[inline(always)]
    pub fn get_charge(&self) -> &i8 {
        &self.charge
    }

    #[inline(always)]
    pub fn get_charge_mut(&mut self) -> &mut i8 {
        &mut self.charge
    }

    #[inline(always)]
    pub fn set_charge(&mut self, charge: i8) {
        self.charge = charge;
    }
}


impl From<Element> for Ion {
    fn from(e: Element) -> Ion {
        Ion::new(e, 0)
    }
}

impl AsRef<Element> for Ion {
    fn as_ref(&self) -> &Element {
        &self.element
    }
}

impl AsMut<Element> for Ion {
    fn as_mut(&mut self) -> &mut Element {
        &mut self.element
    }
}