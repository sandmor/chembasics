use crate::structural::ParserError;

pub fn parse_entry(line: String) -> Result<(), ParserError> {
    if !line.starts_with("M") {
        return Err(ParserError::Syntax);
    }
    Ok(())
}